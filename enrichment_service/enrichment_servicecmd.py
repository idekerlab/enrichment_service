#!/usr/bin/env python

import os
import sys
import argparse
import json

import numpy as np
import pandas
from ndex2.cx2 import RawCX2NetworkFactory, CX2Network
from gprofiler import GProfiler

import enrichment_service


def _parse_arguments(desc, args):
    """
    Parses command line arguments
    :param desc:
    :param args:
    :return:
    """
    help_fm = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_fm)
    parser.add_argument('input',
                        help='Input: data in node table format.')
    parser.add_argument('--mode',
                        choices=['gprofiler', 'iquery'],
                        default='gprofiler',
                        help='Mode. Default: gprofiler.')
    parser.add_argument('--maxpval', type=float, default=0.00000001,
                        help='Max p value')
    parser.add_argument('--minoverlap', default=0.05, type=float,
                        help='Minimum Jaccard to allow for hits')
    parser.add_argument('--omit_intersections', action='store_true',
                        help='If set, do NOT query for gene intersections')
    parser.add_argument('--excludesource', default='HP,MIRNA,TF',
                        help='Comma delimited list of sources to exclude')
    parser.add_argument('--maxgenelistsize', type=int,
                        default=500, help='Maximum number of genes that can'
                                          'be passed in via a query, '
                                          'exceeding this results in '
                                          'error')
    parser.add_argument('--precision', type=int, default=3,
                        help='Number of decimal places to round '
                             'jaccard')
    parser.add_argument('--organism', default='hsapiens',
                        help='Organism to use')
    return parser.parse_args(args)


def run_gprofiler(genes, maxgenelistsize, organism, maxpval, omit_intersections, minoverlap, excludesource, precision,
                  gprofwrapper=GProfiler(user_agent='enrichment-service/' + enrichment_service.__version__,
                                         return_dataframe=True)):
    """
    todo
    """
    if ',' in genes:
        genes = genes.strip(',').strip('\n').split(',')
    else:
        genes = genes.strip(' ').strip('\n').split(' ')
    genelist_size = len(genes)
    if genes is None or genelist_size == 0 or (genelist_size == 1 and len(genes[0].strip()) == 0):
        sys.stderr.write('No genes found in input')
        return None
    if genelist_size > maxgenelistsize:
        sys.stderr.write('Gene list size of ' +
                         str(genelist_size) +
                         ' exceeds max gene list size of ' +
                         str(maxgenelistsize))
        return None

    df_result = gprofwrapper.profile(query=genes, domain_scope="known",
                                     organism=organism,
                                     user_threshold=maxpval,
                                     no_evidences=omit_intersections)

    if not isinstance(df_result, pandas.DataFrame) or df_result.shape[0] == 0:
        return None

    df_result['Jaccard'] = 1.0 / (1.0 / df_result['precision'] +
                                  1.0 / df_result['recall'] - 1)

    # filter out any rows where min overlap is not met
    df_result.drop(df_result[df_result.Jaccard < minoverlap].index,
                   inplace=True)

    # filter out any rows in source exclude list
    if excludesource is not None:
        for a_source in excludesource.split(','):
            df_result.drop(df_result[df_result.source == a_source].index,
                           inplace=True)

    if df_result.shape[0] == 0:
        return None

    # sort by Jaccard and fallback to p_value
    df_result.sort_values(['Jaccard', 'p_value'],
                          ascending=[False, True], inplace=True)

    df_result.reset_index(drop=True, inplace=True)

    df_result = df_result[:1]
    df_result = df_result.round({'Jaccard': precision})

    annotated_members = df_result['intersections'][0]

    theres = {'CD_CommunityName': df_result['name'][0],
              'CD_AnnotatedMembers': ' '.join(annotated_members),
              'CD_AnnotatedMembers_Size': len(annotated_members),
              'CD_AnnotatedMembers_Overlap': len(annotated_members) / len(genes),
              'CD_AnnotatedMembers_Pvalue': df_result['p_value'][0],
              'CD_Labeled': len(df_result['name'][0]) > 0,
              'CD_AnnotatedAlgorithm': 'gProfiler',
              'CD_NonAnnotatedMembers': ' '.join(list(set(genes) - set(annotated_members))),
              'CD_AnnotatedMembers_SourceDB': df_result['source'][0],
              'CD_AnnotatedMembers_SourceTerm': df_result['native'][0]
              }

    return theres


def run_enrichment_with_gprofiler(node_table, theargs):
    results_for_rows = {}
    for node_id, node_val in node_table["rows"].items():
        genes = node_val["CD_MemberList"]
        results_for_rows[node_id] = run_gprofiler(genes, theargs.maxgenelistsize, theargs.organism, theargs.maxpval,
                                                  theargs.omit_intersections, theargs.minoverlap,
                                                  theargs.excludesource, theargs.precision)

    theres = {
        "columns": [{"id": "CD_CommunityName", "type": "string"},
                    {"id": "CD_AnnotatedMembers", "type": "string"},
                    {"id": "CD_AnnotatedMembers_Size", "type": "integer"},
                    {"id": "CD_AnnotatedMembers_Overlap", "type": "string"}, #TODO
                    {"id": "CD_AnnotatedMembers_Pvalue", "type": "string"}, #TODO
                    {"id": "CD_Labeled", "type": "boolean"}, #TODO
                    {"id": "CD_AnnotatedAlgorithm", "type": "string"},
                    {"id": "CD_NonAnnotatedMembers", "type": "string"},
                    {"id": "CD_AnnotatedMembers_SourceDB", "type": "string"},
                    {"id": "CD_AnnotatedMembers_SourceTerm", "type": "string"}],
        "rows": results_for_rows
    }
    return theres


def read_inputfile(inputfile):
    """

    :param inputfile:
    :return:
    """
    with open(inputfile, 'r') as f:
        return json.load(f)


def main(args):
    """
    Main entry point for program

    :param args: command line arguments usually :py:const:`sys.argv`
    :return: 0 for success otherwise failure
    :rtype: int
    """
    desc = """
    TODO
    """

    theargs = _parse_arguments(desc, args[1:])
    try:
        theres = None

        if theargs.mode == 'gprofiler':
            json_input = read_inputfile(theargs.input)
            theres = run_enrichment_with_gprofiler(json_input, theargs)

        if theres is None:
            sys.stderr.write('No results\n')
        else:
            json.dump(theres, sys.stdout, indent=2)
        sys.stdout.flush()
        return 0
    except Exception as e:
        sys.stderr.write('Caught exception: ' + str(e))
        return 2


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))

#!/usr/bin/env python

import sys
import argparse
import json
import time

import pandas
import requests
from gprofiler import GProfiler

import enrichment_service

SOURCES_KEY = 'sources'
RESULTS_KEY = 'results'
DETAILS_KEY = 'details'
SIMILARITY_KEY = 'similarity'


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
    parser.add_argument('--url', default='http://ndexbio.org',
                        help='Endpoint of REST service')
    parser.add_argument('--polling_interval', default=1,
                        type=float, help='Time in seconds to'
                                         'wait between '
                                         'checks on task '
                                         'completion')
    parser.add_argument('--timeout', default=30,
                        type=int, help='Timeout for http '
                                       'requests in seconds')
    parser.add_argument('--retrycount', default=180, type=int,
                        help='Number times to check for completed'
                             'request. Take this value times'
                             'the --polling_interval to determine'
                             'how long this tool will wait'
                             'for a completed result')
    return parser.parse_args(args)


def get_completed_result(resturl, taskid, user_agent,
                         timeout=30):
    """

    :param resultasdict:
    :return:
    """
    res = requests.get(resturl + '/integratedsearch/v1/' + taskid +
                       '',
                       headers={'Content-Type': 'application/json',
                                'User_agent': user_agent},
                       timeout=timeout)
    if res.status_code != 200:
        sys.stderr.write('Received http error: ' +
                         str(res.status_code) + '\n')
        return None
    return res.json()


def wait_for_result(resturl, taskid, user_agent, polling_interval=1,
                    timeout=30,
                    retrycount=180):
    """
    Polls **resturl** with **taskid**
    :param resturl:
    :param taskid:
    :param user_agent:
    :param polling_interval:
    :param timeout:
    :param retrycount:
    :return: True if task completed successfully False otherwise
    :rtype: bool
    """
    counter = 0
    while counter < retrycount:
        try:
            res = requests.get(resturl + '/integratedsearch/v1/' +
                               taskid + '/status',
                               headers={'Content-Type': 'application/json',
                                        'User_agent': user_agent},
                               timeout=timeout)

            if res.status_code == 200:
                jsonres = res.json()
                if jsonres['progress'] == 100:
                    if jsonres['status'] != 'complete':
                        sys.stderr.write('Got error: ' + str(jsonres) + '\n')
                        return False
                    return True
            else:
                sys.stderr.write('Received error : ' +
                                 str(res.status_code) +
                                 ' while polling for completion')
        except requests.exceptions.RequestException as e:
            sys.stderr.write('Received exception waiting for task'
                             'completion: ' + str(e))

        counter += 1
        time.sleep(polling_interval)
    return False


def get_best_result_by_similarity(resultasdict):
    """
    Gets best result by cosine similarity. The service
    currently sorts by pvalue, but this is not returning
    the best result so lets use cosine similarity

    :param resultasdict:
    :type resultasdict: dict
    :return: best result as a dict
    :rtype: dict
    """
    best_similarity = -1.0
    best_result = None
    for cursource in resultasdict[SOURCES_KEY]:
        for curresult in cursource[RESULTS_KEY]:
            if best_similarity == -1.0 or \
                    curresult[DETAILS_KEY][SIMILARITY_KEY] > best_similarity:
                best_result = curresult
                best_similarity = curresult[DETAILS_KEY][SIMILARITY_KEY]
    return best_result


def get_result_in_mapped_term_json(resultasdict, genes):
    """

    :param resultasdict:
    :return:
    """

    if resultasdict is None:
        sys.stderr.write('Results are None\n')
        return None

    if SOURCES_KEY not in resultasdict:
        sys.stderr.write('No sources found in results\n')
        return None

    if resultasdict[SOURCES_KEY] is None:
        sys.stderr.write('Source is None\n')
        return None

    if len(resultasdict[SOURCES_KEY]) <= 0:
        sys.stderr.write('Source is empty\n')
        return None

    if RESULTS_KEY not in resultasdict[SOURCES_KEY][0]:
        sys.stderr.write('Results not in source\n')
        return None

    if resultasdict[SOURCES_KEY][0][RESULTS_KEY] is None:
        sys.stderr.write('First result is None')
        return None

    if len(resultasdict[SOURCES_KEY][0][RESULTS_KEY]) <= 0:
        sys.stderr.write('No result found\n')
        return None

    bestresult = get_best_result_by_similarity(resultasdict)

    colon_loc = bestresult['description'].find(':')
    if colon_loc == -1:
        source = 'NA'
    else:
        source = bestresult['description'][0:colon_loc]

    annotated_members = bestresult['hitGenes']
    name = bestresult['description'][colon_loc + 1:].lstrip()

    theres = {'CD_CommunityName': name,
              'CD_AnnotatedMembers': ' '.join(annotated_members),
              'CD_AnnotatedMembers_Size': len(annotated_members),
              'CD_AnnotatedMembers_Overlap': len(annotated_members) / len(genes),
              'CD_AnnotatedMembers_Pvalue': bestresult['details']['PValue'],
              'CD_Labeled': len(name) > 0,
              'CD_AnnotatedAlgorithm': 'iQuery',
              'CD_NonAnnotatedMembers': ' '.join(list(set(genes) - set(annotated_members))),
              'CD_AnnotatedMembers_SourceDB': source,
              'CD_AnnotatedMembers_SourceTerm': 'NA'
              }

    return theres


def run_iquery(genes, theargs):
    """
    todo

    :param inputfile:
    :param theargs:
    :param gprofwrapper:
    :return:
    """
    genes = genes.strip(',').strip('\n').split(',')
    if genes is None or (len(genes) == 1 and len(genes[0].strip()) == 0):
        sys.stderr.write('No genes found in input')
        return None
    user_agent = 'enrichment-service/' + enrichment_service.__version__
    resturl = theargs.url

    query = {'geneList': genes,
             'sourceList': ['enrichment']}
    res = requests.post(resturl + '/integratedsearch/v1/',
                        json=query, headers={'Content-Type': 'application/json',
                                             'User-Agent': user_agent},
                        timeout=theargs.timeout)
    if res.status_code != 202:
        sys.stderr.write('Got error status from service: ' + str(res.status_code) + ' : ' + res.text + '\n')
        return None

    taskid = res.json()['id']

    if wait_for_result(resturl, taskid, user_agent,
                       timeout=theargs.timeout,
                       retrycount=theargs.retrycount,
                       polling_interval=theargs.polling_interval) is False:
        return None

    resjson = get_completed_result(resturl, taskid, user_agent,
                                   timeout=theargs.timeout)
    return get_result_in_mapped_term_json(resjson, genes)


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


def run_enrichment(node_table, theargs, mode):
    results_for_rows = {}
    for node_id, node_val in node_table["rows"].items():
        genes = node_val["CD_MemberList"]
        if mode == 'gprofiler':
            results_for_rows[node_id] = run_gprofiler(genes, theargs.maxgenelistsize, theargs.organism, theargs.maxpval,
                                                      theargs.omit_intersections, theargs.minoverlap,
                                                      theargs.excludesource, theargs.precision)
        elif mode == 'iquery':
            results_for_rows[node_id] = run_iquery(genes, theargs)
        else:
            raise Exception("Pick either gprofiler or iquery algorithm")

    theres = {
        "columns": [{"id": "CD_CommunityName", "type": "string"},
                    {"id": "CD_AnnotatedMembers", "type": "string"},
                    {"id": "CD_AnnotatedMembers_Size", "type": "integer"},
                    {"id": "CD_AnnotatedMembers_Overlap", "type": "string"},  # TODO
                    {"id": "CD_AnnotatedMembers_Pvalue", "type": "string"},  # TODO
                    {"id": "CD_Labeled", "type": "boolean"},  # TODO
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

        json_input = read_inputfile(theargs.input)
        theres = run_enrichment(json_input, theargs, theargs.mode)

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

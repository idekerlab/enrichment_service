#!/usr/bin/env python

import os
import sys
import argparse
import json
from ndex2.cx2 import RawCX2NetworkFactory, CX2Network

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
                        help='Input: network in CX2 format, node data or edge data.')
    return parser.parse_args(args)

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

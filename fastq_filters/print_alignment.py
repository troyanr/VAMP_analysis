
import os
import re
import argparse

from .filter_generics import print_alignment
from .base_parsers.filter_base_parser import SEQUENCE_RE

def write_allignment(fname, expected_sequence, line_length=200,
                     exclude_at_start=0, exclude_at_end=0, query_offset=0):
    '''
    Print formatted alignments for sequences in `fname`.
    '''

    # check arguments
    if not (query_offset >= 0 and exclude_at_start >= 0):
        raise RuntimeError('Invalid query offset or exclude length')
    if exclude_at_start > query_offset:
        hit_padding = exclude_at_start - query_offset
        query_offset = 0
    else:
        query_offset = query_offset - exclude_at_start
        hit_padding = 0

    with open(fname, 'r') as inF:
            lines = [x.strip() for x in inF.readlines()]

    ofname = '{}.alignment.txt'.format(os.path.splitext(fname)[0])
    with open(ofname, 'w') as outF:
        for line in lines:
            _query = line[exclude_at_start:]
            _query = _query[:(len(_query) - exclude_at_end)]
            _query = '{}{}'.format(' ' * hit_padding, _query)
            print_alignment(_query, expected_sequence,
                            line_length=line_length,
                            query_offset=query_offset, out=outF,
                            query_name='expected', hit_name='observed')


def main():
    parser = argparse.ArgumentParser(description='Print formatted alignment between expected_sequence and sequences in input_files')

    parser.add_argument('--query_offset', default=0, type=int,
                    help='Offset at which expected_sequence should begin in query.')
    parser.add_argument('--exclude_at_start', default=0, type=int,
            help='Length at beginning of alignment to ignore (length of diversity sequence)')
    parser.add_argument('--exclude_at_end', default=0, type=int,
            help='Exclude any bases at the end? (e.g. because many truncated runs)')
    parser.add_argument('--line_length', default=200, type=int,
            help='Number of characters to print in each line. Default is 200.')

    parser.add_argument('expected_sequence', help='Expected sequence. Must match the regex "{}"'.format(SEQUENCE_RE))
    parser.add_argument('input_files', nargs='+', help='Input txt file(s) containing hit sequences.')

    args = parser.parse_args()

    if re.search(SEQUENCE_RE, args.expected_sequence):
        _expected_sequence = args.expected_sequence
    else:
        raise RuntimeError('{} is an invalid sequence!'.format(args.expected_sequence))

    for fname in args.input_files:
        write_allignment(fname, _expected_sequence, query_offset=args.query_offset,
                         exclude_at_start=args.exclude_at_start, exclude_at_end=args.exclude_at_end,
                         line_length=args.line_length)



import sys
import re
import os
import argparse
from multiprocessing import cpu_count
import time

from Bio import SeqIO

from .filter_generics import FilterResult, run_filter
from .base_parsers import FILTER_BASE_PARSER
from .base_parsers.filter_base_parser import OPTIONS as filter_options
from .base_parsers.filter_base_parser import SEQUENCE_RE
from .constants import EXCLUDE_AT_END, EXCLUDE_AT_START, DEFAULT_Q0,\
                      DEFAULT_Q1, DEFAULT_Q2, DEFAULT_QF


def q_score_filter(query, Q_scores, expected_sequence,
                   query_offset=0,
                   exclude_at_start=EXCLUDE_AT_START,
                   exclude_at_end=EXCLUDE_AT_END,
                   Q1=DEFAULT_Q1,
                   Q2=DEFAULT_Q2,
                   QF=DEFAULT_QF,
                   Q0=DEFAULT_Q0,
                   verbose_filter=False):

    '''
    Determine whether query sequence reads pass Q score thresholds.

    Parameters
    ---------
    query: str
        Sequence given by read.
    Q_scores: list
        Quality scores for query. Should be a list of numeric values the same
        length as query.
    expected_sequence: str
        Expected sequence.
    query_offset: int
        Offset at which expected_sequence should begin in query.
    exclude_at_start: int
        Number of bases to exclude from analysis at the beginning of the
        sequence, typically the diversity sequence
    exclude_at_end: int
        Number of bases to exclude from analysis at the end of the sequence,
        often due to decreased quality scores
    Q1: int
        Minimum allowed score for the full sequence
    QF: int
        Number of bases below Q1 allowed per sequence
    Q0: int
        Absolute minimum allowed score for bases allowed by QF
    Q2: int
        Minimum allowed score for bases which do not match expected_sequence.
    verbose_filter: bool
        Show verbose output for each sequence filtered?

    Returns
    -------
    results: FilterResult
    '''

    _Q_scores = Q_scores[query_offset:]
    _query = query[query_offset:]
    begin_index = exclude_at_start
    end_index = len(expected_sequence) - exclude_at_end

    if verbose_filter:
        print_alignment(expected_sequence, _query, line_length=200)
        __query = ''
        __expected_sequence = ''

    # Check whether lengths match up.
    assert(len(_Q_scores) == len(_query))
    if len(_Q_scores) < end_index or len(_query) < end_index:
        return FilterResult(False, 'len', query)

    failed_bases_count = 0
    for i in range(begin_index, end_index):
        if verbose_filter:
            __query += _query[i]
            __expected_sequence += expected_sequence[i]

        if _Q_scores[i] < Q0:
            return FilterResult(False, 'Q0', _query)
        if _Q_scores[i] < Q1:
            failed_bases_count += 1
        if failed_bases_count > QF:
            return FilterResult(False, 'Q1', _query)
        if _query[i] != expected_sequence[i]:  # check if _query and expected base match
            if _Q_scores[i] < Q2:  # If they don't match, apply Q2 filter.
                return FilterResult(False, 'Q2', query)

    if verbose_filter:
        print(__query)
        print(__expected_sequence + '\n')

    return FilterResult(True, None, query)


def q_score_helper(seq_obj, expected_sequence, **kwargs):
    '''
    Receive a Bio.SeqIO.SeqRecord and call q_score_filter with the appropriate
    parameters.

    Parameters
    ----------
    seq_obj: Bio.SeqIO.SeqRecord
        Individual record to
    expected_sequence: str
        Expected read sequence.
    **kwargs:
        Additional arguments to pass to q_score_helper

    Returns
    -------
    results: FilterResult
    '''
    return q_score_filter(str(seq_obj.seq),
                          seq_obj.letter_annotations['phred_quality'],
                          expected_sequence, **kwargs)


def run_qscore_filter(fastq_path, expected_sequence,
                      output_dir_path=None, **kwargs):
    '''
    Read FASTQ files, filter reads by quality scores, and write passing
    sequences to a text file.

    Parameters
    ---------
    fastq_path: str
        File path of fastq file to filter.
    expected_sequence: str
        Expected sequence for the full read.
    output_dir_path: str
        Path to directory to write output files.
        If None, writes to same directory as fastq_path.
    **kwargs:
        Additional arguments passed to run_filter and q_score_filter.

    Returns
    -------
    output_file_name : str
        Path to file where passing sequences have been written.
    '''

    # initialize output file names
    fastq_base = os.path.splitext(os.path.basename(fastq_path))[0]
    output_dir = os.path.abspath(os.path.dirname(fastq_path)) if output_dir_path is None else output_dir_path
    output_file_name = '{}/{}.qScoreFiltered.txt'.format(output_dir, fastq_base)
    failed_file_name = '{}/{}.qScoreFailed.txt'.format(output_dir, fastq_base)

    # start timer
    start_time = time.time()
    sys.stdout.write('Working on {}...\n'.format(fastq_path))
    sys.stdout.write('\tStarted filtering reads at {}\n'.format(time.strftime('%H:%M:%S',
                                                                time.localtime())))

    # Get fastq entry generator
    fastq_sequences = SeqIO.parse(fastq_path, format='fastq')

    # Run filter
    results_count = run_filter(output_file_name, failed_file_name,
                               fastq_sequences, q_score_helper, expected_sequence,
                               verbose=True, use_pool=False, **kwargs)

    # Print summary statistics to stdout
    out = sys.stdout
    total = results_count.total()
    out.write('\n\tQ score filtering results for {}\n'.format(fastq_base))
    out.write('\t{:,} reads processed\n'.format(total))
    out.write('\t{:,} had short Q_score lists ({:.2%})\n'.format(results_count.results_count['len'],
                                                                 results_count.results_count['len'] / total))
    for filter_name in ('Q0', 'Q1', 'Q2'):
        if filter_name in kwargs:
            filter_value = kwargs[filter_name]
        else:
            filter_value = globals()['DEFAULT_{}'.format(filter_name)]
        out.write('\t{:,} failed {} = {} ({:.2%})\n'.format(results_count.results_count[filter_name],
                                                            filter_name,
                                                            filter_value,
                                                            results_count.results_count[filter_name] / total))
    out.write('\tFiltered results written to {}\n'.format(output_file_name))

    # print elapsed time
    end_time = time.time()
    out.write('Total run time = {:.0f} seconds ({:.1f} minutes)\n'.format(end_time - start_time,
                                                                          (end_time - start_time) / 60))


def main():
    parser = argparse.ArgumentParser(description='Apply Q filter to fastq file(s).',
                                     parents=[FILTER_BASE_PARSER])
    args = parser.parse_args()

    if re.search(SEQUENCE_RE, args.expected_sequence):
        _expected_sequence = args.expected_sequence
    else:
        raise RuntimeError('{} is an invalid sequence!'.format(args.expected_sequence))

    for fname in args.input_files:
        run_qscore_filter(fname, _expected_sequence,
                          nThread=args.nThread,
                          output_dir_path=args.output_dir,
                          chunk_size=1000, update_period=1e5,
                          **{x: getattr(args, x) for x in filter_options})

if __name__ == '__main__':
    main()



import sys
from collections import Counter
import multiprocessing as mp
import functools
from abc import ABC, abstractmethod
from math import ceil
from itertools import zip_longest


class FilterResult(object):
    '''
    FilterResult objects store information about the results of a filter
    applied to a sequence read.
    '''

    __slots__ = ['passed', 'fail_reason', 'sequence']

    def __init__(self, passed, fail_reason, sequence):
        self.passed = passed
        self.fail_reason = fail_reason
        self.sequence = sequence


class BaseResultCount(ABC):

    def __init__(self):
        self.results_count = Counter()

    @abstractmethod
    def update(self, results):
        self.results_count.update([r.fail_reason for r in results])

    def total(self):
        ''' Get the total number of results '''
        return sum(self.results_count.values())

    def print_update(self, out):
        '''
        Print the number of failing and passing sequences to out

        Parameters
        ---------
        out: ostream
            Any object with a write method.
        '''

        failed = 0
        passed = 0
        for k, v in self.results_count.items():
            if k is None:
                passed += v
            else:
                failed += v
        out.write('\t{:,} sequences passed, {:,} failed...\n'.format(passed,
                                                                     failed))


class BasicResultCount(BaseResultCount):
    def update(self, results):
        super().update(results)


class dummy_pool():

    ''' dummy context manager functions '''
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, tracebadk):
        return False

    def map(self, f, *args):
        for x in zip(*args):
            yield f(*x)


def _read_chunk(seq_gen, current_seq_index, next_seq_chunk,
                chunk_size=int(1e6)):
    '''
    Read a chunk of sequence elements from a sequence generator.

    Parameters
    ---------
    file_gen: generator
        generator yielding individual sequence elements.
    current_seq_index: int
        Last index processed.
    next_seq_chunk:
        Last index to process in next chunk.
    chunk_size: int
        Number of elements to read.

    Returns
    -------
    done: bool
        Have we reached the end of seq_gen?
    current_seq_index int:
        Updated current_seq_index
    next_seq_chunk: int
        Updated next_seq_chunk
    ret: list
        List of sequence elements.
    '''

    done = False
    ret = list()
    try:
        for _ in range(current_seq_index, next_seq_chunk):
            ret.append(next(seq_gen))
    except StopIteration:
        done = True
    finally:
        current_seq_index = next_seq_chunk
        next_seq_chunk = current_seq_index + chunk_size

    return done, current_seq_index, next_seq_chunk, ret


def run_filter(output_file_name, failed_file_name,
               seq_generator, filter_fxn, expected_sequence,
               chunk_size=1e3, update_period=1e5, nThread=1,
               verbose=True, results_count=None, use_pool=True, **kwargs):
    '''
    Apply filter_fxn over all sequences in seq_generator.

    Parameters
    ---------
    output_file_name: str
        Path to file to write passing sequences.
    failed_file_name: str
        Path to file to write failing sequences.
    seq_generator: generator
        A generator yielding individual sequence elements.
    filter_fxn: function
        A function which acts on individual sequence elements. The function
        should return a FilterResult object.
    verbose: bool, default True
        Produce verbose output?
    chunk_size: int, default 1e3
        Number of sequences to hold in memory before writing to output file.
    update_period: int, default 1e5
        How often should the progress be updated.
    nThread: int, default 1.
        Number of threads to use for parallel processing.
    results_count: A subclass of BaseResultCount, default None
    use_pool: bool, default True
        Use thread pool to process embarrassingly parallel data?
        
    **kwargs
        Additional arguments passed to filter_fxn.

    Returns
    -------
    results_count: A subclass of BaseResultCount
        An object containing a summary of the results.
        If the results_count argument was left as None, returns a
        BaseResultCount object.
    '''

    # Initialize output files
    for fname in (output_file_name, failed_file_name):
        with open(fname, 'w') as outF:
            outF.write('')

    # Initialize control variables
    current_seq_index = 0
    chunk_size = int(chunk_size)
    next_seq_chunk = chunk_size
    if results_count is None:
        results_count = BasicResultCount()
    done = False
    n_updates = 0

    with mp.Pool(processes=nThread) if use_pool else dummy_pool() as pool:
        while not done:
            # First get a chunk of sequences and store them in a temporary list
            done, current_seq_index, next_seq_chunk, sequences_chunk = _read_chunk(seq_generator,
                                                                                   current_seq_index,
                                                                                   next_seq_chunk,
                                                                                   chunk_size=chunk_size)

            # apply filter function using process pool
            results = list(pool.map(functools.partial(filter_fxn,
                                                      expected_sequence=expected_sequence,
                                                      **kwargs),
                                    sequences_chunk))
            results_count.update(results)

            # Write chunk to output files
            with open(output_file_name, 'a') as good_outF,\
                    open(failed_file_name, 'a') as bad_outF:
                for result in results:
                    if result.passed:
                        good_outF.write('{}\n'.format(result.sequence))
                    else:
                        bad_outF.write('{}\n'.format(result.sequence))

            # Count up failed and passed sequences for update
            if verbose:
                update_num = (current_seq_index - (current_seq_index % update_period)) // update_period
                if update_num > n_updates:
                    n_updates += 1
                    results_count.print_update(sys.stdout)

    return results_count


def print_alignment(hit, query, midline=None,
                    query_name='query', hit_name='hit', midline_name='',
                    match_char='|', mismatch_char=' ', truncated_char=' ',
                    query_offset=0,
                    line_length=100, out=sys.stdout):
    '''
    Print a nicely formatted alignment between `hit and `query` to out.
    
    Parameters
    ---------
    hit: str
        Hit sequence.
    query: str
        Query sequence.
    midline: str, default None
        Midline sequence. If None, it will be calculated automatically.
    query_name: str, default 'query'
    hit_name: str, default 'hit'
    midline_name: str, default 'midline'
    match_char: str, default '|'
        Midline char for matching bases.
    mismatch_char: str, default ' '
        Midline char for mismatching bases.
    truncated_char: str, default ' '
        Midline char for truncated hit or query.
    query_offset: int, default 0
        Offset position of query in hit
    line_length: int, default 100
        Number of characters to print per line.
    out: ostream, default sys.stdout
        Stream to write to.
    '''

    def _write_element(out, name, value):
        out.write('{}\t{}\n'.format(name, value))
    
    def _max_length(sequences):
        return max([len(seq) for seq in sequences])
    
    query = '{}{}'.format(' ' * query_offset, query)

    if midline is None:
        midline = ''
        for q, h in zip_longest(query, hit):
            if q is None or h is None or q == ' ' or h == ' ':
                midline += truncated_char
            elif q != h:
                midline += mismatch_char
            else:
                midline += match_char

    max_name_len = _max_length((query_name, hit_name, midline_name))
    length = _max_length((query, hit, midline))
    n_lines = ceil(length / line_length)
    for i in range(n_lines):
        begin = i * line_length
        end = (i + 1) * line_length
        end = end if end < length else length

        for seq, name in zip((query, midline, hit), (query_name, midline_name, hit_name)):
            seq_end = len(seq[:end].replace(' ', ''))
            print_seq = seq[begin:end]
            if print_seq == '':
                print_seq = ' '
            print_seq += ' ' * ((end - begin) - len(print_seq))
            _write_element(out, '{}{}'.format(name, ' ' * (max_name_len - len(name))),
                           '{} {}'.format(print_seq, seq_end))

        out.write('\n')




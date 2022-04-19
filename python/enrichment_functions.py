
import csv
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

def count_sequences(files, exclude_at_start, expected_length, prefix=None):
    '''
    Count unique sequences in file with filtered reads.

    Parameters
    ----------
    fname: list
        Path to sequence files.
    exclude_at_start: int
        Length at beginning of alignment to ignore (length of diversity sequence)
    expected_length: int
        Exclude any bases at the end? (e.g. because many truncated runs)
    prefix: str
        Optional prefix to add to file names before reading

    Returns
    -------
    ret: Counter
        A Counter populated with the counts of unique sequences.
    '''

    ret = Counter()
    for f in files:
        _fname = f if prefix is None else '{}/{}'.format(prefix, f)
        with open(_fname, 'r') as inF:
            for line in inF:
                line = line.strip()[exclude_at_start:]
                line = line[:expected_length]
                if line in ret:
                    ret[line] += 1
                else:
                    ret[line] = 1

    return ret


def parse_sequence(seq, exclude_at_start, exclude_at_end):
    '''
    Remove residues from beginning and end of sequence

    Parameters
    ----------
    exclude_at_start: int
        Length at beginning of sequence to ignore
    exclude_at_end: int
        Length at end of sequence to ignore

    Returns
    -------
    seq: str
        Modified sequence.
    '''
    seq = seq[exclude_at_start:]
    seq = seq[:len(seq) - exclude_at_end]
    return seq


def read_library_sequences(fname, exclude_at_start, exclude_at_end):
    '''
    Read library input sequences.

    Parameters
    ----------
    fname: str
        Path to file.
    exclude_at_start: int
        Length at beginning of sequence to ignore
    exclude_at_end: int
        Length at end of sequence to ignore

    Returns
    -------
    ret: dict
        Set of library sequences.
    '''
    ret = dict()
    with open(fname, 'r') as inF:
        dialect = csv.Sniffer().sniff(inF.read(1024), delimiters=',\t')
        inF.seek(0)
        reader = csv.reader(inF, dialect)

        for line in reader:
            seq = parse_sequence(line[1], exclude_at_start, exclude_at_end)
            ret[seq] = line[0]

        return ret

def get_mutations(base_seq, mod_seq, offset=0):
    '''
    Get the mutations between `base_seq` and `mod_seq`

    Returns
    -------
    ret: str
        Mutation signature.

    Raises
    ------
    AssertionError
        If the lengths of `base_seq` and `mod_seq` are not the same.
    '''
    assert(len(base_seq) == len(mod_seq))
    
    mutation_indecies = ','.join(['{}{}{}'.format(b, i - offset, m) for i, (b, m) in enumerate(zip(base_seq, mod_seq)) if b != m])
    return mutation_indecies


def make_base_matrix_columns(sequences):
    '''
    Get a unique list of bases found in `sequences`.
    '''
    columns = set()
    for x in sequences:
        columns.update(set(x))
    return columns


def populate_logo_matrix(base_seq, sequences, quantifications, columns=None):
    '''
    Populate sequence logo matrix with quantification at each position.

    Parameters
    ----------
    base_seq: str
        Unmodified base sequence.
    sequences: list
        List of library sequences.
    quantification: list
        List of quantification for each library sequence.
    columns: list
        List of unique bases which appear in library.
        If None, make_base_matrix_columns is called to make the column names.

    Returns
    -------
    ret: pd.DataFrame
        Data fram with columns for each base and rows for each base position.

    Raises
    ------
    AssertionError
        If the length of sequences and quantification are not the same.
    '''
    
    if columns is None:
        columns = make_base_matrix_columns(sequences + list(base_seq))
    
    # initialize empty quantification matrix
    ret = pd.DataFrame(np.zeros(shape=(len(base_seq), len(columns))),
                       index=list(range(len(base_seq))),
                       columns=columns)

    assert(len(sequences) == len(quantifications))
    for s, q in zip(sequences, quantifications):
        if not np.isnan(q):
            for i, b in enumerate(s):
                if b != base_seq[i]:
                    ret[b][i] += q

    return ret


class LogoFeature():
    __slots__ = ['begin_index', 'end_index', 'color']

    def __init__(self, begin_index=0, end_index=0, color='b'):
        self.begin_index = int(begin_index)
        self.end_index = int(end_index)
        self.color = color

    def __repr__(self):
        return 'LogoFeture({}, {}, {})'.format(self.begin_index,
                                               self.end_index,
                                               repr(self.color))

    def set_color_rgb(self, rgb):
        if any([x > 1 for x in rgb]):
            self.color = tuple([x / 256 for x in rgb])
        else:
            self.color = tuple(rgb)


def feature_cordinates(base_seq, feature_seq, label_offset=0):
    ''' Get the coordinates of a feature in the parent sequence.'''
    begin = base_seq.find(feature_seq)
    if begin == -1:
        return None, None
    return begin + label_offset, begin + len(feature_seq) - 1 + label_offset


def make_logo_plot(fname, base_seq, dat, width=14, height=4, features=None):
    '''
    Make sequence logo plot.

    Parameters
    ----------
    fname: str
        Path to write logo file.
    base_seq: str
        Sequence to show in plot
    dat: DataFrame
        4xN matrix with numbers to show for each base at each position.
    width: int
        Width in inches.
    height: int
        Height in inches.
    features: list
        List of LogoFeatures to add to plot.
    '''

    x_labels = list()
    base_colors = {'A': 'g', 'C': 'b', 'G': 'darkorange', 'T': 'r'}
    for base, index in zip(base_seq, dat.index.to_list()):
        if index % 5 == 0:
            x_labels.append('{}\n{}'.format(base, index))
        else:
            x_labels.append(base)

    # create logo plot
    logo = logomaker.Logo(dat)
    logo.ax.set_ylabel('$\log_2$ fold enrichment')
    logo.ax.set_xticks(dat.index.to_list())
    logo.ax.set_xticklabels(x_labels)
    logo.style_glyphs_below(flip=False)
    for b, t in zip(base_seq, logo.ax.get_xticklabels()):
        t.set_color(base_colors[b])

    # add feature annotations
    if features is not None:
        for f in features:
            logo.highlight_position_range(pmin=f.begin_index, pmax=f.end_index, color=f.color)

    logo.fig.set_size_inches((width, height))
    logo.fig.savefig(fname, format='pdf')
    plt.close('all')


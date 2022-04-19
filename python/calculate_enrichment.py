
import sys
import argparse
import os
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enrichment_functions import *
from read_params import read_params

def check_files(files, name, prefix=None):
    '''
    Check whether files exist and print results.

    Parameters
    ----------
    files: list
        List of files
    name: str
        Sample name
    prefix: str
        Optional prefix path to add to file paths

    Returns
    -------
    all_found: bool
        True if all files were found.
    '''
    n_files = len(files)
    all_found = n_files > 0
    print('{} {} files:'.format(n_files, name))
    for f in files:
        found = os.path.isfile('{}/{}'.format(data_dir, f))
        fount_s = 'Found' if found else 'Not found'
        print('\t{} -> {}'.format(f, fount_s))
        if not found:
            all_found = False
    return all_found


# Parse command line arguments
parser = argparse.ArgumentParser(description='Calculate enrichment at positions and make sequence logos.')
parser.add_argument('params_file', help='Path to params.yaml')
args = parser.parse_args()

# parse params file
params = read_params(args.params_file)

# define variables from params
selection_names = list(params['file_names']['output'].keys())

# setup external parameters
base_seq = params['base_seq'] # sequence to show in plots
sequence_offset = params['sequence_offset']
expected_len = params['expected_len']
library_exclude_at_start = params['library_exclude_at_start']
library_exclude_at_end = params['library_exclude_at_end']
data_dir = params['path_to_data_dir'] # path to dir with filtered sequence files
read_threshold_in_input = params['read_threshold_in_input'] # Minimum number of reads for a sequence to be used.
library_seq_path = params['path_to_library_sequences']
aggregate_samples = None if 'aggregate_samples' not in params else params['aggregate_samples']
files = dict()
files['input'] = params['file_names']['input']
for sample, sample_files in params['file_names']['output'].items():
    files[sample] = sample_files

# output file names
enriched_sequences_name = params['enriched_sequences_name']
logo_prefix = params['logo_prefix']
matrix_prefix = params['matrix_prefix']
logo_postfix = params['logo_postfix']
matrix_postfix = params['matrix_postfix']
logo_features = params['logo_features']

# check sample names
for name in selection_names:
    if aggregate_samples is not None and name in aggregate_samples:
        print('ERROR: Aggregate name "{}" is also a selection name!')
        sys.exit(-1)
    if name == 'input':
        print('ERROR: Cannot use "input" as a selection name!')
        sys.exit(-1)
if aggregate_samples is not None:
    if 'input' in aggregate_samples:
        print('ERROR: Cannot use "input" as an aggregate name!')
        sys.exit(-1)
    for samples in aggregate_samples.values():
        for sample in samples:
            if sample not in selection_names + ['input']:
                print('ERROR: Requested aggregate sample "{}" is not a known sample!')
                sys.exit(-1)

# check files
all_found = all([check_files(flist, name, data_dir) for name, flist in files.items()])
if not all_found:
    print('ERROR: Not all requested input files exist!')
    sys.exit(-1)
if aggregate_samples is not None:
    print()
    for agg_name, samples in aggregate_samples.items():
        print('Aggregating samples to "{}":'.format(agg_name))
        for sample in samples:
            print('\t{}'.format(sample))

# count sequences
dat = pd.DataFrame()
sums = dict()
for sample_name, lanes in files.items():
    sequences = count_sequences(lanes, sequence_offset, expected_len, prefix=data_dir)
    dat = dat.append(pd.DataFrame(data={'file': [sample_name for _ in range(len(sequences))],
                                        'sequence': list(sequences.keys()),
                                        'count': list(sequences.values())}),
                     ignore_index=True)
    sums[sample_name] = sum([x for x in sequences.values()])

# apply read_threshold_in_input filter
before_len = len(dat.index)
print('\nApplying read_threshold filter:')
dat = dat[dat.apply(lambda x: True if x['file'] != 'input' else x['count'] > read_threshold_in_input, axis=1)]
after_len = len(dat.index)
print('\t{} sequences before, {} after.'.format(before_len, after_len))

# remove sequences which were not found in the input
input_sequences = set(dat[dat['file'] == 'input']['sequence'].to_list())
dat = dat[dat['sequence'].apply(lambda x: x in input_sequences)]

# add fraction of total column
dat['fraction_of_total'] = dat['count'] / dat.groupby(['file'], group_keys=False)['count'].transform(sum)

# add fold enrichment
input_fractions = dict()
for i, r in dat[dat['file'] == 'input'].iterrows():
    assert(r['sequence'] not in input_fractions)
    input_fractions[r['sequence']] = r['fraction_of_total']
dat['fold_enrichment'] = dat.apply(lambda x: x['fraction_of_total'] / input_fractions[str(x['sequence'])] if
        x['sequence'] in input_fractions else np.inf, axis=1)

# add enrichment factor column
dat['enrichment_factor'] = dat['fold_enrichment'] / dat.groupby(['file'], group_keys=False)['fold_enrichment'].transform(max)

# add enrichment rank column
dat['enrichment_rank'] = dat.groupby('file')['fold_enrichment'].rank('dense', ascending=False)

# Get library sequences
lib_seq = read_library_sequences(library_seq_path, library_exclude_at_start, library_exclude_at_end)

# add library rank
dat['in_library'] = dat['sequence'].apply(lambda seq: seq in lib_seq)
dat['library_enrichment_rank'] = dat[dat['in_library']].groupby('file')['fold_enrichment'].rank('dense', ascending=False)

# pivot by sequence
pivot = dat.pivot(index='sequence', columns='file',
                  values=['count', 'fraction_of_total', 'fold_enrichment', 'enrichment_factor', 'enrichment_rank', 'library_enrichment_rank'])
pivot = pivot[pivot.columns[pivot.columns.reindex(['input'] + selection_names, level='file')[1]]]

# Add sequence info columns
pivot.insert(0, 'modifications', pd.Series({x: get_mutations(base_seq, x, offset=3) for x in pivot.index}))
pivot.insert(0, 'library_accession', pd.Series({x: lib_seq[x] if x in lib_seq else np.nan for x in pivot.index}))
pivot.insert(0, 'in_library', pd.Series({x: x in lib_seq for x in pivot.index}))

# add aggregate columns
if aggregate_samples is not None:
    for name, samples in aggregate_samples.items():
        pivot['avg_enrichment', name] = pivot['enrichment_factor'][samples].apply(np.mean, axis=1)
        pivot['stdev', name] = pivot['enrichment_factor'][samples].apply(np.std, axis=1)
        pivot['avg_rank', name] = pivot['avg_enrichment', name].rank(method='dense', ascending=False)
        pivot['avg_library_rank', name] = pivot[pivot['in_library']]['avg_enrichment', name].rank(method='dense', ascending=False)

pivot.to_csv(enriched_sequences_name, sep='\t')

# make seq logo matrix
mat_cols = make_base_matrix_columns(pivot.index.to_list() + list(base_seq))
input_mat = populate_logo_matrix(base_seq, pivot.index.to_list(), pivot['fraction_of_total', 'input'].to_list(), mat_cols)
logo_dats = dict()
agg_f = np.vectorize(lambda n, d: 0 if n == 0 else math.log2(n / d))
for sele in selection_names:
    mat = populate_logo_matrix(base_seq, pivot.index.to_list(), pivot['fraction_of_total', sele].to_list(), mat_cols)
    logo_dats[sele] = pd.DataFrame(agg_f(mat, input_mat), columns=mat_cols)
    logo_dats[sele].index.name = 'position'
    logo_dats[sele].index = list(range(1, expected_len + 1))

# add aggregate logo matrices
stdev_dats = dict()
if aggregate_samples is not None:
    for name, samples in aggregate_samples.items():
        np.mean([mat for name, mat in logo_dats.items() if name in samples], axis = 0)
        logo_dats[name] = pd.DataFrame(np.mean([mat for name, mat in logo_dats.items() if name in samples], axis = 0), columns=mat_cols)
        logo_dats[name].index.name = 'position'
        logo_dats[name].index = list(range(1, expected_len + 1))
        stdev_dats[name] = pd.DataFrame(np.std([mat for name, mat in logo_dats.items() if name in samples], axis = 0), columns=mat_cols)
        stdev_dats[name].index.name = 'position'
        stdev_dats[name].index = list(range(1, expected_len + 1))

# print stdev_dats
for k, v in stdev_dats.items():
    v.to_csv('{}{}.stdev.{}.tsv'.format(matrix_prefix, k, logo_postfix), sep = '\t')

# make sequence logos
for name, mat in logo_dats.items():
    # save logo matrix
    mat.to_csv('{}{}.{}.tsv'.format(matrix_prefix, name, logo_postfix), sep = '\t')

    make_logo_plot('{}{}.{}.pdf'.format(logo_prefix, name, logo_postfix), base_seq, mat)
    make_logo_plot('{}{}.{}.only_positive.pdf'.format(logo_prefix, name, logo_postfix),
                   base_seq, mat.applymap(lambda x: x if x > 0 else 0))

    logo_annotations = list()
    for f, c in logo_features.items():
        cord = feature_cordinates(base_seq, f, 1)
        logoFeature = LogoFeature(cord[0], cord[1])
        if isinstance(c, tuple):
            logoFeature.set_color_rgb(c)
        else:
            logoFeature.color = c
        logo_annotations.append(logoFeature)

    make_logo_plot('{}{}.{}.pdf'.format(logo_prefix, name, logo_postfix), base_seq, mat, features=logo_annotations)
    make_logo_plot('{}{}.{}.only_positive.pdf'.format(logo_prefix, name, logo_postfix),
                   base_seq, mat.applymap(lambda x: x if x > 0 else 0),
                   features=logo_annotations)

plt.close('all')


# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 17:20:15 2019

@author: Rachel Kelemen

This is a brief example piece of code to make a "tile 2" promoter library
"""

import csv
import argparse

def library_generator(promoter_seq, possible_bases,
                      exclude_at_start=0, exclude_at_end=0,
                      double_tile=False):
    """
    A library generator function

    Inputs
    ---
    promoter_seq : str
    possible_bases : list
    exclude_at_start: int
        Number of bases to exclude from analysis at the beginning of the
        sequence, typically the diversity sequence
    exclude_at_end: int
        Number of bases to exclude from analysis at the end of the sequence,
        often due to decreased quality scores
    double_tile : bool

    Returns
    ---
    list_of_sequences : list
        [sequence_name, sequence]
    """

    # Check the promoter sequence for violations
    # Stops the loop if an invalid character (not A, C, G, or T) is found
    for c in promoter_seq:
        if c not in possible_bases:
            print('Error\nlibrary_generator ended because of invalid character in promoter_seq')
            return None

    # Set up a table to hold the resulting sequences
    list_of_sequences = []

    # Add the original sequence
    current_sequence = promoter_seq
    current_name = 'input_seq'
    list_of_sequences.append([current_name, current_sequence])

    sequence_count = 0
    if double_tile:
        # Double tile
        for i in range(exclude_at_start, len(promoter_seq) - exclude_at_end):
            # Same as single tile, but need to consider 2 bases
            # We only want to add sequences where both positions are different, so
            # make two lists of possible bases, one for each position, with the current base
            # removed.
            # Loop through all possible combinations of possible bases (9 total) and add
            # to the list of sequences.
            a, b = promoter_seq[i], promoter_seq[i+1]
            temp_possible_bases_a = possible_bases.copy()
            temp_possible_bases_b = possible_bases.copy()
            temp_possible_bases_a.remove(a)
            temp_possible_bases_b.remove(b)
            for m in temp_possible_bases_a:
                for n in temp_possible_bases_b:
                    current_sequence = promoter_seq[:i] + m + n + promoter_seq[i+2:]
                    current_name = '{}-{}_{}-{}'.format(i+1, m, i+2, n)
                    list_of_sequences.append([current_name, current_sequence])
                    sequence_count += 1
    else:
        # Single tile
        for i in range(exclude_at_start, len(promoter_seq) - exclude_at_end):
            # Identify the current base in the promoter
            # Create a copy of the list of possible bases and remove the current base
            # Generate sequences with each other possible base at that position
            # Add those sequences and their names to the list
            c = promoter_seq[i]
            temp_possible_bases = possible_bases.copy()
            temp_possible_bases.remove(c)
            for n in temp_possible_bases:
                current_sequence = promoter_seq[:i] + n + promoter_seq[i+1:]
                current_name = '{}-{}'.format(i+1, n)
                list_of_sequences.append([current_name, current_sequence])
                sequence_count += 1

    print(str(sequence_count) + ' sequences generated from single tiling')

    return list_of_sequences


def write_to_csv(list_of_sequences, file_path):
    """
    Takes the list of names and sequences from the previous function and writes them
    to a .csv file using the Python csv module

    Inputs
    ---
    list_of_sequences : list
    file_path : str

    Returns
    ---
    None
    """

    # Open the file (this creates the .csv file if it does not already exist)
    # Opening the file in 'w' (write) mode automatically erases any data already present
    # Using the 'with open' statement closes the file when the loop concludes
    # Write each [name, sequence] as one row
    with open(file_path, 'w', newline='') as output_file:
        output_writer = csv.writer(output_file)
        for sequence in list_of_sequences:
            output_writer.writerow(sequence)

def main():
    parser = argparse.ArgumentParser(description='Generate tiled library sequences.')

    parser.add_argument('-b', '--bases', default = 'ACGT',
                        help='Bases to allow as single string. Default is "ACGT".')
    parser.add_argument('-o', '--ofname', default = 'library.csv',
                        help='Name of output file. Default is "library.csv"')
    parser.add_argument('--exclude_at_start', default=0, type=int,
                        help='Length at beginning of alignment to ignore (length of diversity sequence)')
    parser.add_argument('--exclude_at_end', default=0, type=int,
                        help='Exclude any bases at the end? (e.g. because many truncated runs)')
    parser.add_argument('--double', default=False, action='store_true',
                        help='Use double variable bases.')

    parser.add_argument('library_sequence', help='Base sequence used to generate library.')

    args = parser.parse_args()

    # Promoter sequence must be in all caps with T, not U
    possible_bases = [b for b in args.bases]

    list_of_sequences = library_generator(args.library_sequence, possible_bases,
                                          exclude_at_start=args.exclude_at_start,
                                          exclude_at_end=args.exclude_at_end,
                                          double_tile=args.double)

    print('\n{:,} total sequences generated'.format(len(list_of_sequences)))

    write_to_csv(list_of_sequences, args.ofname)
    print('\nSequences written to')
    print(args.ofname)


if __name__ == '__main__':
    main()


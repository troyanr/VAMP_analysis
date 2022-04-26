
import argparse

OPTIONS = ['randomized_length', 'additional_mismatches']

MISMATCH_FILTER_BASE = argparse.ArgumentParser(add_help=False)

MISMATCH_FILTER_BASE.add_argument('--randomized_length', default=2, type=int,
                                  help='Length of randomized region.')

MISMATCH_FILTER_BASE.add_argument('--additional_mismatches', default=0, type=int,
                                  help='Number of additional mismatches to allow.')


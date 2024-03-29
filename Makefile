########## Parameters which will need to be edited for each analysis ##########

# Parameters for mismatch and q score filters
QUERY_OFFSET := --query_offset 34
EXCLUDE_LEN := --exclude_at_start 25 --exclude_at_end 25
FILTER_OPTIONS := $(QUERY_OFFSET) --nThread 1 $(EXCLUDE_LEN)
MISMATCH_FILTER_OPTIONS := --randomized_length 1

# Starting library sequence
EXPECTED_SEQUENCE := CCGCCTAGGAGTACGGTCTCGCTCGAATATTTGCATGTCGCTATGTGTTCTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATCAGACCACAGATCCCCGGAACGAGACCCGATCCGCTCGCAC

# files produced by calculate_enrichment script
FINAL_DATA := tables/enriched_sequences.tsv
 # := tables/enriched_sequences.tsv tables/logo_matrix

###############################################################################

# Get the names of all target filtered fastq sequences
FILTERED_FILES = $(shell ls data/*.fastq|sed 's/\.fastq/.qScoreFiltered.mismatchFiltered.txt/')

# Paths to executable fastq filter scripts
Q_SCORE_FILTER := venv/bin/q_score_filter
MISMATCH_FILTER := venv/bin/mismatch_filter
PRINT_ALLIGNMENT := venv/bin/print_alignment
FASTQ_FILTER_EXE := $(Q_SCORE_FILTER) $(MISMATCH_FILTER) $(PRINT_ALLIGNMENT)

.PHONY: all clean-code clean-data clean-all check-and-reinit-submodules

all: code data

data: $(FILTERED_FILES) data/library_sequences.csv $(FINAL_DATA)

code: venv/touchfile venv $(FASTQ_FILTER_EXE)

# Init fastq_filters git submodules
check-and-reinit-submodules:
	@if git submodule status | egrep -q '^[-]|^[+]' ; then \
		echo "INFO: Need to reinitialize git submodules"; \
		git submodule update --init; \
	fi

venv: venv/touchfile

$(FASTQ_FILTER_EXE): venv
	make check-and-reinit-submodules
	./venv/bin/pip install -e fastq_filters

# Init python virtual environment
venv/touchfile: requirements.txt
	test -d venv || python3 -m venv venv
	./venv/bin/pip install -Ur requirements.txt
	touch venv/touchfile

clean-code:
	rm -rf venv

clean-data:
	rm -fv data/*.txt
	rm -fv data/library_sequences.csv
	rm -fv tables/* fig/*

clean-all: clean-code clean-data

# optional target to print alignment files
ALIGNMENT_FILES = $(shell ls data/*|egrep '(Failed|Filtered).txt$$'|sed 's/\.txt/.alignment.txt/')
alignment: $(ALIGNMENT_FILES)
data/%.alignment.txt: $(PRINT_ALLIGNMENT) data/%.txt
	$< $(QUERY_OFFSET) $(EXPECTED_SEQUENCE) $(word 2,$^)

# run q score filter
data/%.qScoreFiltered.txt: $(Q_SCORE_FILTER) data/%.fastq
	$< $(FILTER_OPTIONS) $(EXPECTED_SEQUENCE) $(word 2,$^)

# run mismatch filter
data/%.qScoreFiltered.mismatchFiltered.txt: $(MISMATCH_FILTER) data/%.qScoreFiltered.txt
	$< $(FILTER_OPTIONS) $(MISMATCH_FILTER_OPTIONS) $(EXPECTED_SEQUENCE) $(word 2,$^)

# build library sequences
data/library_sequences.csv: python/generate_tiled_library.py
	./venv/bin/python3 python/generate_tiled_library.py $(EXCLUDE_LEN) -o data/library_sequences.csv $(EXPECTED_SEQUENCE)

# aggregate data and make sequence logos
$(FINAL_DATA): python/calculate_enrichment.py python/enrichment_functions.py python/read_params.py params.yaml $(FILTERED_FILES)
	mkdir -p fig tables
	./venv/bin/python python/calculate_enrichment.py params.yaml

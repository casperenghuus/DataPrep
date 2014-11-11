import sys
import os
import re

import NGS_CNE
import RNA_contigs

DEFAULT_REGEX = 's_G1_L001_R([1-2])_([0-9]+).fastq.([0-9]+).gz'
# CWD = os.getcwd()
CWD = '/scratch/cne/ecre/'

target = re.compile('202.fixed.fa')
for root, dirs, files in os.walk(CWD, followlinks=True):
    for file_string in files:
        file_match = target.search(file_string)
        if file_match:
            REF_SEQ_FILE = os.path.join(root, file_string)


def run_seqprepRNA(file1, file2,
    default_regex=DEFAULT_REGEX,
    CWD=CWD,
    ref_seq_file=REF_SEQ_FILE):
    '''
    Requires NGS_CNE.py to be in the same folder as this file.
    1) A file pair is entered as forward and reverse files using bash.
    it should be entered as a relative path.
    2) A file of all the reference sequences is made
    3) Run SeqPrep on the file pair. Adapters are used.
    4) Uses bash to extract sequences from the SeqPrep results
    5) Trims the merged sequences
    6) Compares the merged sequence to the trimmed reference sequences.
    if there is an exact match, the frequency is counted, and the
    sequence is added to a new file.
    7) Provide some basic 'statistics' on the output
    '''

    # Get bin and file-numbers.
    file_regex = re.compile(default_regex)
    bin_number = file_regex.search(file1).group(3)
    file_number = file_regex.search(file1).group(2)

    print '\n### Initiating SeqPrep #'
    print '### Paired files:', file1, '<>', file2

    # Get reference sequences from 202.fixed.fa. Make sure the file DO NOT
    # contain (trace spacer) in the fasta header (or any other white-space
    # character for that matter)!!!
    has_header = True
    dictionary = False
    NGS_CNE.load_reference_fasta(ref_seq_file, has_header, dictionary)

    list_to_be_tested = NGS_CNE.load_fq_files()

    # Sort to bins to test for coherence:
    NGS_CNE.sort_to_bins(DEFAULT_REGEX, list_to_be_tested)

    # Run SeqPrep
    RNA_contigs.initiate_seqprep(file_number, bin_number, file1, file2)

    # Run statistics
    RNA_contigs.write_stats(file_number, bin_number, file1)

if __name__ == "__main__":

    # Run SeqPrep taking forward file and
    # reverse file as the two required inputs.
    run_seqprepRNA(*sys.argv[1:])

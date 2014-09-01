#import sys
import os
import re # may be omitted

import NGS_CNE

# Set paths/files
# CWD = os.getcwd()
CWD = '/mnt/sda1/My-Documents/Dropbox/casper_ecre/'
file_output1 = os.path.join(
    CWD, 'results/sequences.allresults.txt')
file_output2 = os.path.join(
    CWD, 'results/stats.allresults.txt')
# REF_SEQ_FILE = os.path.join(CWD, 'fa/202.fixed.fa')

target = re.compile('202.fixed.fa')
for root, dirs, files in os.walk(CWD, followlinks=True):
    for file_string in files:
        file_match = target.search(file_string)
        if file_match:
            REF_SEQ_FILE = os.path.join(CWD, file_string)

# Set parameters
has_header = True
dictionary = True
file_regex = '(s.*).results.([0-9]+).([0-9]+).*'

if __name__ == "__main__":

    # Get reference sequences
    name_dict = NGS_CNE.load_reference_fasta(
        REF_SEQ_FILE, has_header, dictionary)

    # Creates a list of stat results to be merged
    sequence_list, stats_list = NGS_CNE.generate_file_list(file_regex)

    sequence_list = NGS_CNE.sort_to_bins(file_regex, sequence_list)
    stats_list = NGS_CNE.sort_to_bins(file_regex, stats_list)

    # Test that the two lists have equal length. Should always be True
    assert len(sequence_list) == len(stats_list)

    # Go through each bin, and merge the files from these. 
    NGS_CNE.merge_all(name_dict, sequence_list, stats_list, file_output1, file_output2)

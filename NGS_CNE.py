import re
import os.path
import subprocess
# import sys
import operator

from Bio import SeqIO
from decimal import *
from tabulate import tabulate

# Loads default paths
# CWD = os.getcwd()
# FASTQ_PATH = os.path.join(CWD, 'fq/')
# FASTA_PATH = os.path.join(CWD, 'fa/')
# OUTPUT_PATH = os.path.join(CWD, 'gz/')
# RESULTS_PATH = os.path.join(CWD, 'results/')
# REF_SEQ_FILE = os.path.join(CWD, 'fa/202.fixed.fa')

FQ_DIR = '/scratch/cne/ecre/fq/202_hiseq/' #CHANGE!
CWD = '/scratch/cne/ecre/'
# CWD = os.getcwd()
# FQ_DIR = os.path.join(CWD, 'RNA_test/')

DEFAULT_REGEX = '^s_G1_L001_R([1-2])_([0-9]+).fastq.([0-9]+).gz'
DEFAULT_RESTRICTION_SITE_1 = 'CATATG'
DEFAULT_RESTRICTION_SITE_2 = 'GGCGCGCC'
DEFAULT_ADAPTER2 = 'CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT'
DEFAULT_ADAPTER1 = 'GGCGCGCCATGACTAAGCTTTTCATTGTCATGC'
# READ_TRIM_REGEX = '^(.*{restriction_site1})?(.*?)({restriction_site2}.*)?$'
READ_TRIM_REGEX = '(.*{restriction_site1})?(.*)({restriction_site2}){1,2}|({restriction_site2})?'

def load_fq_files(fq_dir=FQ_DIR, target_regex=DEFAULT_REGEX):
    '''
    Load files from folder and stores them in a list of tuples. the
    tuples contain the forward and reverse file from the same bin and
    with the same file number.
    '''
    # Required lists
    file_pairs = []
    forward_list = []
    reverse_list = []

    # Regex target to grab input files.
    target = re.compile(target_regex)

    # Searches through folders and files in the fq_dir. Grabs files
    # matching the regex
    for root, dirs, files in os.walk(fq_dir, followlinks=True):
        for file_string in files:
            file_match = target.search(file_string)
            if file_match:

                # Whether it is a forward (=1) or reverse (=2) file
                direction = int(file_match.group(1))
                if direction == 1:
                    forward_list.append(os.path.join(root, file_string))
                if direction == 2:
                    reverse_list.append(os.path.join(root, file_string))

    # Grab for each file in the forward list, it goes through each file in the
    # reverse list. If it is a long list, it may be worthwhile to use k.pop()
    # to speed up the search.
    for i in forward_list:
        i = os.path.split(i)[1]
        for k in reverse_list:
            k = os.path.split(k)[1]
            # Grabs file from the same bin
            if int(target.search(i).group(2)) == int(target.search(k).group(2)):
                # Grabs file with the same 'sample number'
                if int(
                    target.search(i).group(3)) == int(
                        target.search(k).group(3)):
                    # Appends the pairs as a tuple, (forward, reverse), to a
                    # list of file-pair-tuples.
                    file_pairs.append(
                        (os.path.join(fq_dir, i), os.path.join(fq_dir, k)))

    return file_pairs

def load_reference_fasta(
    ref_seq_file,
    has_header,
    dictionary_instead_of_file,
    read_trim_regex=READ_TRIM_REGEX,
    restriction_site1=DEFAULT_RESTRICTION_SITE_1,
    restriction_site2=DEFAULT_RESTRICTION_SITE_2):
    '''
    Loads reference fasta file and prepare the revere complement of the
    sequence.
    If header = True, the header will be included in the output file.
    Otherwise it will be omitted.
    If dictionary_instead_of_file is True, a file will NOT be made. Instead,
    a dictionary is prepared, which holds the sequence as key and the
    description as value.
    '''

    # Regex used to trim sequence. Case-insensitive
    regex_target = re.compile(read_trim_regex.format(
        restriction_site1=restriction_site1,
        restriction_site2=restriction_site2),
        flags=re.IGNORECASE)

    # A file is made rather than a dictionary over the results. Can
    # be with or without header. Notice the difference in suffix.
    if dictionary_instead_of_file is False:

        # Defines the name of the output file
        if has_header is True:
            trimmed_ref_seq_file = os.path.join(CWD, 'fa/202.trimmed.fixed.fa')
        else:
            trimmed_ref_seq_file = os.path.join(CWD, 'fa/202.trimmed.fixed.seq')

        # Write to output file
        with open(trimmed_ref_seq_file, 'w') as fh:
            for seq_record in SeqIO.parse(ref_seq_file, 'fasta'):
                # Write all sequences in capital letters
                # and get the reverse complement
                seq_record.seq = str(
                    seq_record.seq.upper().reverse_complement())
                # Trim the sequence
                trim = regex_target.search(seq_record.seq)

                # If the sequence has been trimmed, write trimmed sequence
                # to file. Include header if set to True
                if trim is not None:
                    if has_header is True:
                        fh.write('>'+seq_record.id+'\n')
                    fh.write(trim.group(2)+'\n')
                # Else, write the un-trimmed sequence to file. Include
                # header of True
                else:
                    if has_header is True:
                        fh.write('>'+seq_record.id+'\n')
                    fh.write(seq_record.seq+'\n')

        return trimmed_ref_seq_file

    # A dictionary is made containing the sequence as key and
    # fasta header as value.
    if dictionary_instead_of_file is True:

        name_dict = {}

        for seq_record in SeqIO.parse(ref_seq_file, 'fasta'):
            # Write all sequences in capital letters
            # and get the reverse complement
            seq_record.seq = str(seq_record.seq.upper().reverse_complement())
            # Trim the sequence
            trim = regex_target.search(seq_record.seq)

            # Write trimmed sequence to dictionary with header as value
            if trim is not None:
                name_dict[trim.group(2)] = seq_record.id
            # Write un-trimmed sequence to dictionary with header as value
            else:
                name_dict[seq_record.seq] = seq_record.id

        return name_dict

def sort_to_bins(regex, file_list, number_position=2, bin_position=3):
    '''Sort files to bins and performs a coherence check to see if there
    are files missing. I.e. tests if the file numbers are in sequential order
    '''
    # Make an epmty list to be updated with a list of bins
    bin_list = []
    # Regex target
    target = re.compile(regex)

    # Assigns the maximum number of bins
    for i in range(1, len(file_list)+1):
        bin_i = []
        # Goes through each file in the file list,
        # and appends it to a bin-number
        for f in file_list:
            if int(target.search(str(f)).group(int(bin_position))) == int(i):
                bin_i.append(str(f))
        bin_list.append(bin_i)
        if bin_i == []:
            raise IOError('Bin {bin_number} is missing'.format(bin_number=i))

        # Test to see if all bins have files in sequential order.
        number_list = []
        for entry in bin_i:
            file_number = int(target.search(entry).group(int(number_position)))
            # If the file number is not in the list, it is appended to it. Only
            # unique numbers are kept (i.e. forward and reverse files would
            # generate duplicates)
            if file_number not in number_list:
                number_list.append(file_number)

        # Test to see if numbers are in sequential order
        coherence_test = sorted(number_list)

        # This is only True if the list is orderer in steps of +1
        # (it generates a duplicate sorted list used for a True/False statement)
        results = coherence_test == range(
            coherence_test[0], coherence_test[-1]+1)
        # Raises error if the list is not in sequential order
        if results is False:
            raise IOError(
                'Missing files in bin {bin_number}'.format(bin_number=i))

        if len(file_list) == sum([len(bins) for bins in bin_list]):
            break

    assert len(file_list) == sum([len(bin) for bin in bin_list]), (
        'File list length and resulting bin list size do not match')

        # # Tests if all files are used and break out of loop if soin
        # length = 0
        # for k in bin_list:
        #     length += len(k)
        #     if length == len(file_list):
        #         # Breaks all for loops!
        #         # Alternative to return:
        #         # http://stackoverflow.com/questions/653509/breaking-out-of-nested-loops
    return bin_list

def initiate_seqprep(file_number, bin_number, file1, file2,
                adapter1=DEFAULT_ADAPTER1, adapter2=DEFAULT_ADAPTER2):
    '''
    Run SeqPrep using forward and reverse dictionaries as input. Output files
    are located in th same directory as the Handling_NGS_files.py. Output
    files are named OUTPUT_<input file name>
    '''

    # Set name prefixes for outputs
    output_f_prefix = '/scratch/cne/ecre/counts/202_hiseq/output.'+os.path.split(
        file1)[1] # CHANGE!
    output_r_prefix = '/scratch/cne/ecre/counts/202_hiseq/output.'+os.path.split(
        file2)[1] # CHANGE!

    # Run SeqPrep. Refer to manual for the different parameters.
    # Adapters are used.
    subprocess.check_call([
        'SeqPrep',
        '-f', file1,
        '-r', file2,
        '-1', output_f_prefix,
        '-2', output_r_prefix,
        '-3', output_f_prefix[:-3]+'.disc.fq.gz',
        '-4', output_r_prefix[:-3]+'.disc.fq.gz',
        '-s', '/scratch/cne/ecre/counts/202_hiseq/output.'+file_number+'.'+bin_number+'.M.fq.gz',
        '-A', adapter1, '-B', adapter2,
        '-X', '1', '-g', '-L', '5']) #CHANGE

def trim_fq(merged_file_path, restriction_site1=DEFAULT_RESTRICTION_SITE_1,
            restriction_site2=DEFAULT_RESTRICTION_SITE_2, read_trim_regex=READ_TRIM_REGEX):
    '''
    Trim sequences located in same directory as Handling_NGS_files.py (default)
    (the path is defined outside this function!).
    Yields the trimmed sequence (creates a generator)
    '''
    # Regular expression for forward and reverse targets
    regex_target = re.compile(read_trim_regex.format(
        restriction_site1=restriction_site1,
        restriction_site2=restriction_site2),
        flags=re.IGNORECASE)

    zcat_handle = subprocess.Popen(['zcat', merged_file_path],
                                    stdout=subprocess.PIPE)

    while True:

        # grab 4 lines for each record, break from while loop if
        # the iterator is empty at a new header
        try:
            header1 = zcat_handle.stdout.next()
        except StopIteration:
            break

        # However, if iterator is empty inside a record, then we
        # want to to raise the error.
        try:
            seq = zcat_handle.stdout.next()
            header2 = zcat_handle.stdout.next()
            quality = zcat_handle.stdout.next()
        except StopIteration:
            raise IOError('Incomplete record at end of file {}'.format(
                merged_file_path))

        # first attempt to match forward, then try reverse
        trim = regex_target.search(seq)
        if trim is not None:
            yield trim.groups()
        else:
            raise ValueError('Invalid input seuence, trim regex failed!')

    # Used so stdout read the entire line and not just single characters
    zcat_handle.wait()

def grep_merged_read(
    trimming_generator,
    sequence_output,
    stats_output,
    file1,
    restriction_site1=DEFAULT_RESTRICTION_SITE_1,
    restriction_site2=DEFAULT_RESTRICTION_SITE_2):
    '''
    Takes    the trimming_generator as input and matches the trimmed sequences to
    the 202.trimmed.fixed.fa file.
    Output this to a count file, were each unique, matched sequence is sorted
    by its frequency.

    Useful parameters and how to calculate them. Not all are included in output.
    -Number of DNA seq total = sequence_count
    -Number of DNA seq discarded = (zcat disc.fq.gz | wc -l) / 4
    -Number of matched DNA seqs = matched_seq_count
    -Number of unique sequences with no match = absent_seq_tracker
    -Number of unique sequences with a match = present_seq_tracker
    -Percent of total sequences used = matched_seq_count / wrong_seq_count * 100
    -Sequence coverage = present_seq_tracker / len(id_seq_dict.keys())

    '''
    # Iterators
    total_reads = 0
    left_trims = 0
    right_trims = 0
    trimmed_reads = 0
    not_trimmed_reads = 0
    correct_counts = 0

    # Sets significant numbers for calculations (using decimal package)
    getcontext().prec = 3

    # File path is predefined to take the file that has no headers
    trimmed_ref_seq_file = os.path.join(CWD, 'fa/202.trimmed.fixed.seq')

    # Command to grep a perfect matches between the merged SeqPrep file
    # and the file containing the reference sequences. The end result is
    # sorted so only unique sequences and their frequency appears. They
    # are sorted in descending order of frequency
    command = ' | '.join([
        'grep -ixFf {reference_library}',
        'sort',
        'uniq -c',
        'sort -nr > {output_file}']).format(
            reference_library=trimmed_ref_seq_file,
            output_file=sequence_output)

    # Runs the command in bash using the trimmed sequences from the
    # merged SeqPrep file as input (stdin)
    perfect_match = subprocess.Popen(
        command,
        shell=True,
        stdin=subprocess.PIPE)

    # Writes the matched sequences to file and counts the total number
    # of reads as well as the total number of untrimmed sequences and
    # left + right trims.
    for trim_left, read, trim_right in trimming_generator:
        perfect_match.stdin.write(read+'\n')
        total_reads += 1
        if restriction_site1 in read:
            not_trimmed_reads += 1
        elif restriction_site2 in read:
            not_trimmed_reads += 1
        else:
            trimmed_reads += 1
            if trim_left is not None:
                left_trims += 1
            if trim_right is not None:
                right_trims += 1

    # This closes the child process when we are finished sending trimmed seqs
    perfect_match.communicate()

    with open(sequence_output) as seq_handle:
        for line in seq_handle:
            frequency = line.rsplit()[0]
            correct_counts += int(frequency)

    # Writes stats to file. Should be made a separate function...
    with open(stats_output, 'w') as fh:

        # CALCULATIONS/VARIABLES REQUIRED BEFORE WRITING TO FILE:
        # Total number of unmerged sequences:
        # Assumes equal number of lines in forward and reverse file
        unmerged_files = os.path.join(
            CWD, 'counts/202_hiseq/output.'+os.path.split(file1)[1]) #CHANGE!
        command = 'zcat {input_file} | wc -l'.format(
            input_file=unmerged_files)
        unmerged_seqs = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True)
        unmerged_seqs_number = int(unmerged_seqs.communicate()[0])/4

        # Total number of discarded sequences:
        # Assumes equal number of lines in forward and reverse file
        disc_seq_file = unmerged_files[:-3]+'.disc.fq.gz'
        command = 'zcat {input_file} | wc -l'.format(input_file=disc_seq_file)
        disc_seqs = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True)
        disc_seqs_number = int(disc_seqs.communicate()[0])/4

        # Total number of sequences in input file:
        total_seq_num = total_reads + unmerged_seqs_number + disc_seqs_number

        # Percentage of used sequence:
        sequences_used = (
            Decimal(trimmed_reads) +
            Decimal(not_trimmed_reads))/Decimal(total_seq_num)*100

        # Number of unique sequences:
        command = 'cat {input_file} | wc -l'.format(input_file=sequence_output)
        correct_seqs = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True)
        correct_seqs_number = int(correct_seqs.communicate()[0])

        # Percentage of sequences used compared to
        # the total number of sequences given as input:
        percent_correct = Decimal(correct_counts)/Decimal(total_reads)*100

        # Number of reference sequences:
        command = 'cat {input_file} | wc -l'.format(
            input_file=trimmed_ref_seq_file)
        ref_seqs = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        ref_seqs_number = int(ref_seqs.communicate()[0])

        # WRITTEN TO FILE:
        # Total number of sequences:
        fh.write('Total number of sequences given as input:\t'+str(
            total_seq_num)+'\n')

        # Total number of merged sequences
        fh.write('Total merged sequences:\t'+str(total_reads)+'\n')

        # Total number of trimmed sequences:
        fh.write(
            'Total trimmed and merged sequences:\t'+str(trimmed_reads)+'\n')

        # Total number of untrimmed sequences:
        fh.write(
            'Total untrimmed and merged sequences:\t'+str(
                not_trimmed_reads)+'\n')

        fh.write('Percent of total sequences used:\t'+str(sequences_used)+'%\n')

        # Number of left/right trims:
        fh.write('Total left trims ({restriction_site1}):\t'.format(
            restriction_site1=restriction_site1)+str(left_trims)+'\n')
        fh.write('Total right trims ({restriction_site2}):\t'.format(
            restriction_site2=restriction_site2)+str(right_trims)+'\n')

        fh.write('Total unmerged sequences:\t'+str(unmerged_seqs_number)+'\n')
        fh.write('Total discarded sequences:\t'+str(disc_seqs_number)+'\n')
        fh.write(
            'Total unique correct sequences:\t'+str(correct_seqs_number)+'\n')

        # Number of correct and unique sequences:
        fh.write('Total count of correct sequences:\t'+str(correct_counts)+'\n')
        fh. write('Percentage of correct sequences used:\t'+str(
            percent_correct)+'%\n')
        fh.write('Total reference sequences:\t'+str(ref_seqs_number)+'\n')

def generate_file_list(get_files_regex, CWD=CWD):
    '''Generates a list of files based on the regex given. Used to
    prepare a list of files when merging the result summaries'''

    sequence_list = []
    stats_list = []

    # Specifies the folder where statistics may be found
    CWD = os.path.join(CWD, 'counts/202_hiseq') #CHANGE!
    
    target = re.compile(get_files_regex)

    # Prepares a list of files
    for root, dirs, files in os.walk(CWD, followlinks=True):
        for file_string in files:
            file_match = target.findall(file_string)
            if file_match:
                file_path = os.path.join(root, file_string)
                file_match = target.search(file_string)
                if 'sequence' in file_match.group(1):
                    sequence_list.append(file_path)
                if 'stats' in file_match.group(1):
                    stats_list.append(file_path)
    assert len(sequence_list) == len(stats_list)
    return sequence_list, stats_list

def merge_all(name_dict,
    sequence_list,
    stats_list,
    file_output1,
    file_output2,
    restriction_site1=DEFAULT_RESTRICTION_SITE_1,
    restriction_site2=DEFAULT_RESTRICTION_SITE_2):
    '''Merge summaries given for each individual run through SeqPrep'''

    # Sets significant numbers for calculations (using decimal package)
    getcontext().prec = 3

    # Iterator.
    # Reuse of correct_counts variable (it is reset here)
    correct_counts = 0
    total_unique_sequences = 0
    unique_sequences = {}

    # Used to assert that there are an equal number of unique keys and values
    assert len(set(name_dict.values())) == len(set(name_dict.values()))

    # Header for file:
    header = ['Name', 'Total']

    # Make a dictionary which will hold
    # the sum of frequencies across all bins
    total_frequencies = name_dict.fromkeys(name_dict, 0)
    total_nf_dict = {}

    # List of bin dictionaries
    bin_seq_storage = {}

    # Set starting bin number
    i = 0
    # Go through each bin
    for file_table in sequence_list:
        # Append bin number to header
        header.append('Bin.{bin_number},'.format(bin_number=i+1))

        # Make an empty dictionary for bin i to
        # store the name and their frequency
        bin_seq_storage[i] = {}

        # Make another empty dictionary to store unique sequeces for each bin
        unique_sequences[i] = {}

        # For each bin, reset the sequence:frequency dictionary
        bin_frequencies = name_dict.fromkeys(name_dict, 0)

        # Go through every file in a bin
        for file_entry in file_table:
            with open(file_entry) as entry_handle:
                for line in entry_handle:
                    # Split each line to grab the frequency and the sequence
                    count, name = line.rsplit()

                    # Update frequency for bin i
                    bin_frequencies[name] += int(count)

                    # Update frequency for total bins
                    total_frequencies[name] += int(count)

                    # Count the sum of all frequncies across all bins
                    correct_counts += int(count)

                    unique_sequences[i][name] = 1

        # Update the dictonary holding each bin. Each bin consists of a
        # second dictionary of the type name:frequency
        for key_name in total_frequencies.keys():
            bin_seq_storage[i][
            name_dict[key_name]] = bin_frequencies[key_name]

        # Go to next dictionary
        i += 1

    # Update the dictionary over total
    # frequencies to the type name:frequency
    for key_name in total_frequencies.keys():
        total_nf_dict[name_dict[key_name]] = total_frequencies[key_name]

    # Sort the dictionary holding the total number of frequencies
    # in descending order
    sorted_total_freq = sorted(
        total_nf_dict.iteritems(),
        key=operator.itemgetter(1),
        reverse=True)

    # Sets list to be updated with key and frequencies
    bin_seq_values = []

    # Go through each key (sequence name)
    for k in range(0, len(sorted_total_freq)):
        # Sequence name
        key_handle = sorted_total_freq[k][0]
        # Row for each sequence name and its values across bins
        row = []
        # Add sequence name to first column
        row.append(key_handle)
        # Add the total frequency to the second column
        row.append(str(sorted_total_freq[k][1]))
        if sorted_total_freq[k][1] > 0:
            total_unique_sequences += 1
        # Go through each bin in sequential order and append its
        # frequency to the next available column
        for bn in range(0, len(bin_seq_storage.keys())):
            if key_handle in bin_seq_storage[bn].keys():
                row.append(str(bin_seq_storage[bn][key_handle]))

        # Append the finished row to bin_values
        bin_seq_values.append(row)

    # Keys to be used in stats_dict
    left_trim = 'Total left trims ({restriction_site1}):'.format(
        restriction_site1=restriction_site1)
    right_trim = 'Total right trims ({restriction_site2}):'.format(
        restriction_site2=restriction_site2)

    # Dictionary to be updated. The complex keys matches the lines written
    # in the statistics files generated for each merged SeqPrep files
    stats_dict = {
    'Total number of sequences given as input:' : 0,
    'Total merged sequences:' : 0,
    'Total trimmed and merged sequences:' : 0,
    'Total untrimmed and merged sequences:' : 0,
    'Percent of total sequences used:' : 0,
    left_trim : 0,
    right_trim : 0,
    'Total unmerged sequences:' : 0,
    'Total discarded sequences:' : 0,
    'Total unique correct sequences:' : 0,
    'Total count of correct sequences:' : 0,
    'Percentage of correct sequences used:' : 0,
    'Total reference sequences:' : 0}

    # List of bin dictionaries
    bin_stat_storage = {}

    # Set starting bin number
    i = 0

    # Go through each bin
    for file_table in stats_list:
        # Make an empty dictionary for bin i to
        # store the name and their frequency
        bin_stat_storage[i] = {}.fromkeys(stats_dict, 0)
        # Go through every file in a bin
        for file_entry in file_table:
            with open(file_entry) as entry_handle:
                for line in entry_handle:
                    # Grab the text
                    text = line.rsplit('\t')[0]
                    # Grab the number/value
                    value = line.rsplit('\t')[1]
                    # Update dictionary by adding the value to
                    # the current value of the key (=text)
                    if '%' not in value:
                        bin_stat_storage[i][text] += int(value)
            # Retrieve the total number of unique sequences. A dictionary is
            # used, as it removes dupliates
            bin_stat_storage[i][
                'Total unique correct sequences:'] = len(
                unique_sequences[i].keys())

            # Fixes the percentage values as they cannot simply be added:
            bin_stat_storage[i]['Percent of total sequences used:'] = str(
                Decimal(bin_stat_storage[i]['Total merged sequences:'])/Decimal(
                bin_stat_storage[i]['Total number of sequences given as input:'
                ])*100)+'%'
            # More percentages to be fixed
            bin_stat_storage[i][
                'Percentage of correct sequences used:'] = str(Decimal(
                bin_stat_storage[i]['Total count of correct sequences:'
                ])/Decimal(bin_stat_storage[i]['Total merged sequences:'])*100
                )+'%'
            # Predefined number of reference sequences
            bin_stat_storage[i][
                'Total reference sequences:'] = 12653

        # Go to next dictionary
        i += 1

    # List of keynames present in dictionary. Preferred over those in
    # stats_dict, as these are ordered.
    keynames = (
            'Total number of sequences given as input:',
            'Total merged sequences:',
            'Total trimmed and merged sequences:',
            'Total untrimmed and merged sequences:',
            'Percent of total sequences used:',
            left_trim,
            right_trim,
            'Total unmerged sequences:',
            'Total discarded sequences:',
            'Total unique correct sequences:',
            'Total count of correct sequences:',
            'Percentage of correct sequences used:',
            'Total reference sequences:')

    # Make new dictionary to hold the total statistics. Has all keys
    # set to 0
    total_stats_dict = {}.fromkeys(stats_dict, 0)

    # Iterate over each bin
    for i in range(0, len(stats_list)):
        # Iterate over each key
        for kn in keynames:
            if '%' not in str(bin_stat_storage[i][kn]):
                # Take the sum over all bins
                total_stats_dict[kn] += bin_stat_storage[i][kn]
    # Fixes unique sequences using counter
    total_stats_dict[
        'Total unique correct sequences:'] = total_unique_sequences
    # Reuses the predefineed number of reference sequences
    total_stats_dict['Total reference sequences:'] = bin_stat_storage[i][
        'Total reference sequences:']
    # Fixes the percentage values as they cannot simply be added:
    total_stats_dict['Percent of total sequences used:'] = str(
        Decimal(total_stats_dict['Total merged sequences:'])/Decimal(
        total_stats_dict['Total number of sequences given as input:'
        ])*100)+'%'
    total_stats_dict['Percentage of correct sequences used:'] = str(
        Decimal(total_stats_dict['Total count of correct sequences:'
        ])/Decimal(total_stats_dict['Total merged sequences:'])*100
        )+'%'

    # Make empty list to hold the result for each key over every bin
    bin_stat_values = []
    # Go through each key
    for k in range(0, len(keynames)):
        key_handle = keynames[k]
        row = [key_handle+',']
        row.append(str(total_stats_dict[key_handle])+',')
        # Go through each bin
        for bn in range(0, len(bin_stat_storage.keys())):
            row.append(str(bin_stat_storage[bn][key_handle])+',')
        bin_stat_values.append(row)

    # Sequences
    with open(file_output1, 'w') as fh:
        fh.write(tabulate(bin_seq_values, header, tablefmt='plain'))
    # Statistics
    with open(file_output2, 'w') as fh:
        fh.write(tabulate(bin_stat_values, header, tablefmt='plain'))
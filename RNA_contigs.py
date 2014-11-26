import re
import os.path
import subprocess
import operator

import NGS_CNE

from Bio import SeqIO
from decimal import *
from tabulate import tabulate

# Default settings for ssh-server
FQ_DIR = '/scratch/cne/ecre/fq/202_hsdna/' #CHANGE!
CWD = '/scratch/cne/ecre/'
COUNTS_DIR = os.path.join(CWD, 'counts/202_hsdna/')

# Default paths for local computer
# CWD = os.getcwd() # CHANGE!
# FQ_DIR = CWD # CHANGE!

DEFAULT_REGEX = '^s_G1_L001_R([1-2])_([0-9]+).fastq.([0-9]+).gz'
OUTPUT_REGEX = '^output.[0-9]+.([0-9])+.M.fq.gz'
DEFAULT_RESTRICTION_SITE_1 = 'CATATG'
DEFAULT_RESTRICTION_SITE_2 = 'GGCGCGCC'
DEFAULT_ADAPTER_RNA_A = 'CGCCATGACTAAGCTTTTCATTGTC'
DEFAULT_ADAPTER_RNA_B = 'CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT'
DEFAULT_ADAPTER_DNA_A = 'GGCGCGCCATGACTAAGCTTTTCATTGTCATGC' # GGCG removed from 5'-end
DEFAULT_ADAPTER_DNA_B = 'CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT'
READ_TRIM_REGEX = '^(.*{restriction_site1})?(.*?)({restriction_site2}.*)?$'

BIN1_OUTPUT = os.path.join(COUNTS_DIR, 'bin1_counts.seq')
BIN2_OUTPUT = os.path.join(COUNTS_DIR, 'bin2_counts.seq')
ALL_BINS = os.path.join(COUNTS_DIR, 'all_bins.fa')
REF_FASTA = os.path.join(CWD, 'fa/202.trimmed.fixed.fa')
UNMAPPED = os.path.join(COUNTS_DIR, 'unmapped.fa')
BOWTIE_OUT = os.path.join(COUNTS_DIR, 'bowtie_output.csv') 


def initiate_seqprep(file_number, bin_number, file1, file2,
        counts_dir=COUNTS_DIR, adapter_A=DEFAULT_ADAPTER_DNA_A,
        adapter_B=DEFAULT_ADAPTER_DNA_B):
    '''
    Run SeqPrep using forward and reverse dictionaries as input. Output files
    are located in the same directory as the Handling_NGS_files.py. Output
    files are named OUTPUT_<input file name>
    '''

    # # Set name prefixes for outputs
    output_f_prefix = counts_dir+'output.'+os.path.split(file1)[1] # CHANGE!
    output_r_prefix = counts_dir+'output.'+os.path.split(file2)[1] # CHANGE!

    # output_f_prefix = os.path.join(
    #     CWD, 'RNA_test/Output/output.'+os.path.split(file1)[1]) # CHANGE!
    # output_r_prefix = os.path.join(
    #     CWD, 'RNA_test/Output/output.'+os.path.split(file2)[1]) # CHANGE!

    # Run SeqPrep. Refer to manual for the different parameters.
    # Adapters are used.
    subprocess.check_call([
        '/opt/SeqPrep.dbg/SeqPrep',
        '-f', file1,
        '-r', file2,
        '-1', output_f_prefix,
        '-2', output_r_prefix,
        '-3', output_f_prefix[:-3]+'.disc.fq.gz',
        '-4', output_r_prefix[:-3]+'.disc.fq.gz',
        '-s', counts_dir+'output.'+file_number+'.'+bin_number+'.M.fq.gz',
        '-A', adapter_A, '-B', adapter_B,
        '-X', '1', '-g', '-L', '5']) #CHANGE

def seq_counts(output_regex=OUTPUT_REGEX,
    bin1_output=BIN1_OUTPUT,
    bin2_output=BIN2_OUTPUT,
    counts_dir=COUNTS_DIR):
    '''
    1. Make a list of files based on the bin they were sorted to.
    2. Go through one bin at a time
    2.1. zcat multiple files at once = all lines across all files in a bin
    are merged.
    2.2. grep sequences in a way which trims off the adapter.
    2.3. Write to file.
    2.4 Go to next bin
    '''
    # First, make a list of all output files, sorted by bin
    bin_files = [[], []]

    target = re.compile(output_regex)
    for root, dirs, files in os.walk(counts_dir, followlinks=True):
        for f in files:
            file_match = target.search(f)
            if file_match:
                if file_match.group(1) == '1': # DBG = group()
                    bin_files[0].append(os.path.join(counts_dir, f))
                elif file_match.group(1) == '2':
                    bin_files[1].append(os.path.join(counts_dir, f))

    for bin_num in range(0, len(bin_files)):
        # Toggle output file
        if bin_num == 0:
            print '\n>>> BIN1'
            output = bin1_output
        if bin_num == 1:
            print '\n>>> BIN2'
            output = bin2_output

        # Uses zcat on all files for a specific bin
        zcat = 'zcat '+' '.join(bin_files[bin_num])
        # Regex used to grep sequence while manually trimming off adapter
        regex = '^[NATCG]+(?=([NATCG]{2}CGCCATGACTAAGCTTTTCATTGTC))|^[NATCG]+$'
        # This command creates the BinX file
        cmd = "{z} | grep -E '{r}' | sort | uniq -c > {o}".format(
            z=zcat, r=regex, o=output)
        bin_counts = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        bin_counts.communicate()

def merge_bins(ftype, bin1=BIN1_OUTPUT, bin2=BIN2_OUTPUT, all_bins=ALL_BINS):
    '''Merges all sequences for all reads across all bins (uses SeqPrep output)
    A dictionary is made with:
    Key = sequence
    Value[0] = ID#
    Value[1] = Total counts
    Value[2] = Bin1 counts
    Value[3] = Bin2 counts
    '''

    # Get unique sequence while ignoring number. A strange way to do it though
    seq_regex = '(?=[ATCGN])[ATCGN]+'
    # Gets the counts for each unique sequence
    count_regex = '[0-9]+'
    seq = re.compile(seq_regex)
    count = re.compile(count_regex)
    counts_dict = {}

    # Sets necessary REs for trimming DNA
    re1 = re.compile('(?<=CATATG)(.+)', flags=re.IGNORECASE)
    re2 = re.compile('(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)
    re3 = re.compile('(?<=CATATG)(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)

    with open(bin1) as b1:
        print '\n>>>MAKING BIN1'
        # Sets ID = 1 for first entry
        i = 1
        # For each line in Bin1
        for l in b1:
            sequence = str(seq.findall(l)[0])
            counts = int(count.findall(l)[0])

            # THE FOLLOWING CODE GENERATES AN ERROR THROUGH .groups()
            # IF THE REGEX DOESN'T MATCH. IT THEN TRIES THE NEXT ONE 
            # UNTIL THE SEQUENCE IS FULLY TRIMMED.

            # Check if both restriction sites are there.
            # If multiple GGCGCGCC, take only the last
            trimmed = sequence
            try:
                hit1 = re3.search(trimmed)
                hit1.groups()
                trimmed = str(hit1.group(1))

            except:
                # Check if only CATATG is present
                try:
                    hit2 = re1.search(trimmed)
                    hit2.groups()
                    trimmed = str(hit2.group(1))

                except:
                    # Check if only GGCGCGCC is present
                    # If multiple GGCGCGCC, take only the last
                    try:
                        hit3 = re2.search(trimmed)
                        hit3.groups()
                        trimmed = str(hit3.group(1))

                    # No restriction sites present
                    except:
                        trimmed = str(trimmed)

            # Checks if the input type is RNA. Trims off last 2 bp
            if ftype == 'RNA':
                if len(trimmed) > 2:
                    # Trims first two bases of the RNA
                    trimmed = trimmed[:-2]

            # If the sequene is present in the dictionary, add the Bin1
            # and Total counts to the correct spot.
            try:
                counts_dict[trimmed][1] += counts
                counts_dict[trimmed][2] += counts
            # An error is raised if the trimmed sequece is absent.
            # Adds a new entry
            except:
                counts_dict[trimmed] = [
                    str(i), counts, counts, 0]

            i += 1

    with open(bin2) as b2:
        print '\n>>>MAKING BIN2'
        for l in b2:
            sequence = str(seq.findall(l)[0])
            counts = int(count.findall(l)[0])
            
            # THE FOLLOWING CODE GENERATES AN ERROR THROUGH .groups()
            # IF THE REGEX DOESN'T MATCH. IT THEN TRIES THE NEXT ONE 
            # UNTIL THE SEQUENCE IS FULLY TRIMMED.

            # Check if both restriction sites are there.
            # If multiple GGCGCGCC, take only the last
            trimmed = sequence
            try:
                hit1 = re3.search(trimmed)
                hit1.groups()
                trimmed = str(hit1.group(1))

            except:
                # Check if only CATATG is present
                try:
                    hit2 = re1.search(trimmed)
                    hit2.groups()
                    trimmed = str(hit2.group(1))

                except:
                    # Check if only GGCGCGCC is present
                    # If multiple GGCGCGCC, take only the last
                    try:
                        hit3 = re2.search(trimmed)
                        hit3.groups()
                        trimmed = str(hit3.group(1))

                    # No restriction sites present
                    except:
                        trimmed = str(trimmed)

            if ftype == 'RNA':
                if len(trimmed) > 2:
                    # Trims first two bases of the RNA
                    trimmed = trimmed[:-2]

            # A KeyError is raised if the value isn't there.
            # Try adds to existing values while except is a quick and 'dirty'
            # way to add a new, unique entry to the dictionary
            try:
                counts_dict[trimmed][1] += counts
                counts_dict[trimmed][3] += counts
            # Adds a new entry if the trimmed sequence isn't there already.
            except:
                counts_dict[trimmed] = [
                    str(i), counts, 0, counts]
                i += 1

    # WRITE FILE WHICH CAN BE CHECKED FOR ERRORS
    with open(all_bins, 'w') as ab:
        print '\n>>>WRITING FILE!'
        for k in counts_dict:
            # [0] = ID
            # [1] = Total counts
            # [2] = Bin 1 counts
            # [3] = Bin 2 counts
            # k = sequence

            # Write to file
            if len(k) > 4:
                header = '\t'.join(
                    ['>'+str(counts_dict[k][0]),
                    str(counts_dict[k][1]),
                    str(counts_dict[k][2]),
                    str(counts_dict[k][3])+'\n'])
                ab.write(header+str(k)+'\n')

def run_bowtie(all_bins=ALL_BINS, ref_fasta=REF_FASTA,
    unmapped=UNMAPPED, bowtie_out=BOWTIE_OUT,
    bowtie_cmd=None, min_read_count=5, max_read_length=200):
    '''Run Bowtie. Dictionary given as input
    '''

    if not bowtie_cmd:
        # remove the newlines from all headers and change them to tabs
        # then feed the two-column file to bowtie, then filter on minimum
        # read count and maximum read length through perl

        # Makes an output file where sequences matching to multiple
        # promoter have been removed (the reads are too short).
        bowtie_cmd = '''
            bowtie -v 3 -l 10 -k 114 -p 16 \
                --norc --best --strata --suppress 2,6 \
                --un {u} -f {ref}  \
                <(perl -pe 's/([^NATGC])\n/$1\t/' {ab} \
                | perl -ne '@l = split; ($l[1] > {min}
                    && length($l[4]) < {max}
                    && (s/\t([NATGC])/\n$1/ && print));')\
                | awk '{{if($8 <=1) print}}' > {bo}'''.format(
                    u=unmapped, ref=ref_fasta, ab=all_bins,
                    min=min_read_count, max=max_read_length,
                    bo=bowtie_out)

    p = subprocess.Popen(bowtie_cmd,
        stdout=subprocess.PIPE,
        shell=True,
        executable='/bin/bash')
    p.wait()

def bowtie_to_file(all_results, bowtie_out=BOWTIE_OUT, ref_fasta=REF_FASTA):
    '''
    Takes the bowtie-output and writes it to sequeces.allresults.xNA.txt.
    '''

    # Dictionary to hold the Promoter--RBS
    bo_dict = {}
    # Header to print to file
    file_header = ['Name', 'Total', 'Bin.1', 'Bin.2']

    # Makes dictionary with key = Promoter--RBS, value = sequence
    # For each sequence in the reference fasta
    for seq_record in SeqIO.parse(ref_fasta, 'fasta'):
        name = str(seq_record.id)
        bo_dict[name] = [name, 0,0,0]

    # Update dictionary with the counts
    with open(bowtie_out) as bo:
        for l in bo:
            l = l.split()

            # Gets the RBS name
            RBS = l[4].split('--')[1]
            # Get numbr of alternative alignments
            alt_aligns = int(l[7])
            # Gets left-based offset
            l_offset = int(l[5])

            # If the left offset is not zero, the squence is not accepted
            if l_offset > 0:
                continue

            # Get only the best alignment (get read of the alignments
            # which wrongly aligns to the other 113 promoters)
            # BBa_J61100 is never aligned perfectly!
            if alt_aligns == 0:
                # Update Total
                bo_dict[str(l[4])][1] += int(l[1])
                # Update Bin 1
                bo_dict[str(l[4])][2] += int(l[2])
                # Update Bin 2
                bo_dict[str(l[4])][3] += int(l[3])

            # Get the Anderson and BBa_J61100 RBSs.
            # They map to 
            if alt_aligns == 1 and RBS == 'BBa_J61100':
                    # Update Total
                    bo_dict[str(l[4])][1] += int(l[1])
                    # Update Bin 1
                    bo_dict[str(l[4])][2] += int(l[2])
                    # Update Bin 2
                    bo_dict[str(l[4])][3] += int(l[3])

    # Sort dictionary based on total counts in descending order
    sorted_dict = sorted(
        bo_dict.values(),
        key=operator.itemgetter(1),
        reverse=True)

    # Write to file
    with open(all_results, 'w') as ar:
        ar.write(tabulate(sorted_dict, file_header, tablefmt='plain'))

def write_stats(file_number, bin_number, file1, counts_dir=COUNTS_DIR):
    '''
    Write stats for each individual seqprep output
    Useful parameters and how to calculate them. Not all are included in output.
    -Number of DNA seq total = sequence_count
    -Number of DNA seq discarded = (zcat disc.fq.gz | wc -l) / 4
    -Number of matched DNA seqs = matched_seq_count
    -Number of unique sequences with no match = absent_seq_tracker
    -Number of unique sequences with a match = present_seq_tracker
    -Percent of total sequences used = matched_seq_count / wrong_seq_count * 100
    -Sequence coverage = present_seq_tracker / len(id_seq_dict.keys())
    '''

    # Table to write after appending numbers to it
    table = [
        ['Total number of sequences given as input:'],
        ['Total merged sequences:'],
        ['Percent merged:'],
        ['Total unmerged sequences:'],
        ['Total discarded sequences:'],
        ['Total reference sequences:']]

    # File paths
    merged_file = os.path.join(
        counts_dir+'output.'+file_number+'.'+bin_number+'.M.fq.gz')
    unmerged_files = os.path.join(
        counts_dir+'output.'+os.path.split(file1)[1])
    disc_seq_file = unmerged_files[:-3]+'.disc.fq.gz'
    stats_output = os.path.join(
        counts_dir+'stats.results.'+file_number+'.'+bin_number+'.M.txt')

    # Sets significant numbers for calculations (using decimal package)
    getcontext().prec = 3

    # Writes stats to file. Should be made a separate function...
    # ... But it's not
    with open(stats_output, 'w') as so:

        # CALCULATIONS/VARIABLES REQUIRED BEFORE WRITING TO FILE:
        # Total reads
        command = 'zcat {i} | wc -l'.format(i=file1)
        fh = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True)
        total_reads = int(fh.communicate()[0])/4
        table[0].append(str(total_reads))

        # Total merged reads
        command = 'zcat {i} | wc -l'.format(i=merged_file)
        fh = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True)
        merged_reads = int(fh.communicate()[0])/4
        table[1].append(str(merged_reads))

        # Percentage merged:
        table[2].append(str(Decimal(merged_reads)/Decimal(total_reads)*100))

        # Total unmerged reads:
        command = 'zcat {i} | wc -l'.format(i=unmerged_files)
        fh = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True)
        table[3].append(str(int(fh.communicate()[0])/4))

        # Total number of discarded sequences:
        command = 'zcat {i} | wc -l'.format(i=disc_seq_file)
        fh = subprocess.Popen(
            command, stdout=subprocess.PIPE, shell=True)
        table[4].append(str(int(fh.communicate()[0])/4))

        # Number of reference sequences:
        table[5].append('12653')

        # WRITE TO FILE
        so.write(tabulate(table, tablefmt='plain'))

def final_stats(all_stats, all_seqs, counts_dir=COUNTS_DIR):
    '''Merge stats across stat-files and get bowtie stats.
    Number of discarded reads cannot be obtained without
    changing the function running bowtie
    '''

    # Sets significant numbers for calculations (using decimal package)
    getcontext().prec = 3

    target = re.compile('stats.results.([0-9]+).([0-9]+).*')
    stat_files = []

    # Prepares a list of files
    for root, dirs, files in os.walk(counts_dir, followlinks=True):
        for file_string in files:
            file_match = target.findall(file_string)
            if file_match:
                file_path = os.path.join(root, file_string)
                stat_files.append(file_path)

    # Table to be updated with values. Used for writing file
    table = [
        ['Total number of sequences given as input:', 0],
        ['Total merged sequences:', 0],
        ['Percent merged:', 0],
        ['Total unmerged sequences:', 0],
        ['Total discarded sequences:', 0],
        ['Total reference sequences:', 0],
        ['Total aligned sequences', 0],
        ['Total sequences discarded', 0],
        ['Percent of total reads used', 0]]

    with open(all_stats, 'w') as out:
        for sl in stat_files:
            with open(sl) as f:
                i = 0
                for l in f:
                    table[i][1] += float(l.split(':')[1].lstrip().rstrip())
                    i += 1

        # Update some values which are 'unintentionally' changed
        table[2][1] = Decimal(int(table[1][1]))/Decimal(int(table[0][1]))*100
        table[5][1] = '12653'

        # Total aligned sequences:
        cmd = "awk '{{s+=$2}} END {{print s}}' {i}".format(i=all_seqs)
        fh = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        table[6][1] = str(fh.communicate()[0])

        # Total sequences discarded:
        table[7][1] = str(int(table[1][1] - int(table[6][1])))

        # Percent of total reads used:
        table[8][1] = str(
            Decimal(table[6][1])/Decimal(str(table[0][1]))*100)

        # Make everything a string
        for i in range(0, 8):
            table[i][1] = str(table[i][1])

        # Write to file:
        out.write(tabulate(table, tablefmt='plain'))

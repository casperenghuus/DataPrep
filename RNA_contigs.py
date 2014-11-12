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
OUTPUT_REGEX = '^output.[0-9]+.([0-9]).M.fq.gz'
DEFAULT_RESTRICTION_SITE_1 = 'CATATG'
DEFAULT_RESTRICTION_SITE_2 = 'GGCGCGCC'
DEFAULT_ADAPTER1 = 'CGCCATGACTAAGCTTTTCATTGTC'
READ_TRIM_REGEX = '^(.*{restriction_site1})?(.*?)({restriction_site2}.*)?$'

BIN1_OUTPUT = os.path.join(COUNTS_DIR, 'bin1_counts.seq') # CHANGE!
BIN2_OUTPUT = os.path.join(COUNTS_DIR, 'bin2_counts.seq') # CHANGE!
ALL_BINS = os.path.join(COUNTS_DIR, 'all_bins.fa') # CHANGE!
REF_FASTA = os.path.join(CWD, 'fa/202.trimmed.fixed.fa') # CHANGE!
UNMAPPED = os.path.join(COUNTS_DIR, 'unmapped.fa') # CHANGE!
BOWTIE_OUT = os.path.join(COUNTS_DIR, 'bowtie_output.csv') # CHANGE!


def initiate_seqprep(file_number, bin_number, file1, file2, 
        counts_dir=COUNTS_DIR, adapter1=DEFAULT_ADAPTER1):
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
        'SeqPrep',
        '-f', file1,
        '-r', file2,
        '-1', output_f_prefix,
        '-2', output_r_prefix,
        '-3', output_f_prefix[:-3]+'.disc.fq.gz',
        '-4', output_r_prefix[:-3]+'.disc.fq.gz',
        '-s', counts_dir+'output.'+file_number+'.'+bin_number+'.M.fq.gz',
        '-A', adapter1,
        '-X', '1.0', '-g', '-L', '5']) #CHANGE

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
    bins = [[],[]]

    target = re.compile(output_regex)
    for root, files in os.walk(counts_dir, followlinks=True):
        for file_string in files:
            file_match = target.search(file_string)
            if file_match:
                if file_match.group(1) == '1':
                    bins[0].append(os.path.join(root, file_string))
                elif file_match.group(1) == '2':
                    bins[1].append(os.path.join(root, file_string))

    bin_num = 1
    for b in bins:
        # Toggle output file
        if bin_num == 1:
            print '\n>>> BIN1'
            output = bin1_output
        if bin_num == 2:
            print '\n>>> BIN2'
            output = bin2_output

        zcat = 'zcat '+' '.join(b)
        regex = '^[NATCG]+(?=([NATCG]{2}CGCCATGACTAAGCTTTTCATTGTC))|^[NATCG]+$'
        cmd = "{z} | grep -E '{r}' | sort | uniq -c > {o}".format(
            z=zcat, r=regex, o=output)

        bin_counts = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        # Next bin
        bin_num = 2
        bin_counts.communicate()

def merge_bins(bin1=BIN1_OUTPUT, bin2=BIN2_OUTPUT, all_bins=ALL_BINS):
    '''Merges all sequences for all reads across all bins (uses SeqPrep output)
    A dictionary is made with:
    Key = sequence
    Value[0] = ID#
    Value[1] = Total counts
    Value[2] = Bin1 counts
    Value[3] = Bin2 counts
    '''

    seq_regex = '(?=[ATCGN])[ATCGN]+'
    count_regex = '[0-9]+'
    seq = re.compile(seq_regex)
    count = re.compile(count_regex)
    counts_dict = {}

    with open(bin1) as b1:
        print '\n>>>MAKING BIN1'
        i = 1
        for l in b1:
            sequence = str(seq.findall(l)[0])
            counts = int(count.findall(l)[0])

            counts_dict[sequence] = [
                str(i), counts, str(counts), 0]
            i += 1

    with open(bin2) as b2:
        print '\n>>>MAKING BIN2'
        for l in b2:
            sequence = str(seq.findall(l)[0])
            counts = int(count.findall(l)[0])
            # A KeyError is raised if the value isn't there.
            # Try adds to existing values while except is a quick and 'dirty'
            # way to add a new, unique entry to the dictionary
            try:
                counts_dict[sequence][1] += counts
                counts_dict[sequence][3] = str(counts)
            except:
                counts_dict[sequence] = [
                    str(i), counts, 0, str(counts)]
                i += 1

    # WRITE FILE WHICH CAN BE CHECKED FOR ERRORS
    with open(all_bins, 'w') as ab:
        print '\n>>>WRITING FILE!'
        for k in counts_dict:
            header = '\t'.join(
                ['>'+str(counts_dict[k][0]),
                str(counts_dict[k][1]),
                str(counts_dict[k][2]),
                str(counts_dict[k][3])+'\n'])
            ab.write(header+str(k)+'\n')

def run_bowtie(all_bins=ALL_BINS, ref_fasta=REF_FASTA,
    unmapped=UNMAPPED, bowtie_out=BOWTIE_OUT,
    bowtie_cmd=None, min_read_count=1, max_read_length=200):
    '''Run Bowtie. Dictionary given as input
    '''

    if not bowtie_cmd:
        # remove the newlines from all headers and change them to tabs
        # then feed the two-column file to bowtie, then filter on minimum
        # read count and maximum read length through perl

        # CORES ARE SET TO 8! Change to -p 16 for server!!
        # min_read_count HAS BEEN CHANGED! CHANGE BACK TO 5 ON SERVER!!
        # Makes an output file where sequences matching to multiple
        # promoter have been removed (the reads are too short).
        bowtie_cmd = '''
            bowtie -v 3 -l 10 -k 114 -p 16 \
                --norc --best --strata --suppress 2,6 \
                --un {u} -f {ref}  \
                <(perl -pe 's/([^NATGC])\n/$1\t/' {ab} \
                | perl -ne '@l = split; ($l[1] > {min}
                    && length($l[4]) < {max}
                    && (s/\t([NATGC])/\n$1/ && print));') \
                | cut -f -7 | uniq -uf 6 > {bo}'''.format(
                    u=unmapped, ref=ref_fasta, ab=all_bins,
                    min=min_read_count, max=max_read_length,
                    bo=bowtie_out)

    p = subprocess.Popen(bowtie_cmd,
        stdout= subprocess.PIPE,
        shell= True,
        executable='/bin/bash')
    p.wait()

def bowtie_to_file(all_results, bowtie_out=BOWTIE_OUT, ref_fasta=REF_FASTA):
    bo_dict = {}
    file_header = ['Name', 'Total', 'Bin.1', 'Bin.2']

    # Makes dictionary with key = RBS, value = sequence
    for seq_record in SeqIO.parse(ref_fasta, 'fasta'):
        name = str(seq_record.id)
        bo_dict[name] = [name, 0,0,0]

    # Update dictionary with the counts
    with open(bowtie_out) as bo:
        for l in bo:
            l = l.split()
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

    table = [
        ['Total number of sequences given as input:'],
        ['Total merged sequences:'],
        ['Percent merged:'],
        ['Total unmerged sequences:'],
        ['Total discarded sequences:'],
        ['Total reference sequences:']]

    merged_file = os.path.join(
        counts_dir+'output.'+file_number+'.'+bin_number+'.M.fq.gz') # CHANGE!
    unmerged_files = os.path.join(
        counts_dir+'output.'+os.path.split(file1)[1]) #CHANGE!
    disc_seq_file = unmerged_files[:-3]+'.disc.fq.gz'
    
    stats_output = os.path.join(
        counts_dir+'stats.results.'+file_number+'.'+bin_number+'.M.txt') # CHANGE!

    # Sets significant numbers for calculations (using decimal package)
    getcontext().prec = 3

    # Writes stats to file. Should be made a separate function...
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

        # Write to file:
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

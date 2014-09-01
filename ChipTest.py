# NEEDS TO BE DONE:
# A test to examine if every sequence has a primer pair associated with it!

import sys
import os
import re
import subprocess

from Bio import SeqIO

CWD = os.getcwd()
promoters = os.path.join(CWD, 'promoterfile.fasta')
f_primers = os.path.join(CWD, 'skpp15-forward.fasta')
r_primers = os.path.join(CWD, 'skpp15-reverse.fasta')

# res_enz = {
#     bsai : {
#         bsai_cut1_f : 'GGTCTN',
#         bsai_cut2_f : 'CCAGAGNNNNN',
#         bsai_cut1_r : 'NNNNNGAGACC',
#         bsai_cut2_r : 'NCTCTGG'},
#     sapi : {
#         sapi_cut1_f : 'GCTCTTCN',
#         sapi_cut2_f : 'CGAGAAGNNNN',
#         sapi_cut1_r : 'NCTTCTCG',
#         sapi_cut2_r : 'NNNNGAAGAGC'},
#     noti : {noti_f : 'GCGGCCGC', noti_r : 'CGCCGGCG'},
#     ndei : {ndei_f : 'CATATG', ndei_r : 'GTATAC'}}

res_enz = {
    'bsai' : {
        'bsai_cut1_f' : 'GGTCT[ATCG]',
        'bsai_cut2_f' : 'CCAGAG[ATCG]{5}',
        'bsai_cut1_r' : '[ATCG]{5}GAGACC',
        'bsai_cut2_r' : '[ATCG]CTCTGG'},
    'sapi' : {
        'sapi_cut1_f' : '(GCTCTTC.{1})',
        'sapi_cut2_f' : '(CGAGAAG.{4})',
        'sapi_cut1_r' : '(.{1}CTTCTCG)',
        'sapi_cut2_r' : '(.{4}GAAGAGC)'},
    'noti' : 'GCGGCCGC|CGCCGGC',
    'ndei' : {
            'ndei_cut_f' : 'CATATG',
            'ndei_cut_r' : 'CGTATG',
            'ndei_regex' : 'CATATG.+CGTATG'}}

# Table containing sequences with faulty lengths
length_err_header = ['Description']
length_errors = []

res_err_header = ['Description', 'BsaI', 'SapI', 'NotI', 'NdeI']
res_errors = []

### FUNCTIONS:
def length_test(name, sequence, filetype):
    
    # Test for primer length
    if filetype == 'Primer':
        if len(sequence) != 15:
            text = '{fasta_header} is not the designated length'.format(
                fasta_header=name)
            return text    

    # Test for promoter length
    if filetype == 'Promoter':
        if len(sequence) != 230:
            text = '{fasta_header} is not the designated length'.format(
                fasta_header=name)
            return text

def res_site_test(name, sequence, res_enz):
    seq_results = []

    # Test for BsaI:
    res_sites = 0
    for rs in res_enz['bsai']:
        target = re.compile(rs)
        res_sites += len(target.findall(sequence))
    if res_sites == 2:
        bsai = None
    else:
        bsai = 'Error: Wrong number of restriction sites'
    seq_results.append(bsai)

    # Test for SapI restriction site. Must be absent
    for rs in res_enz['sapi']:
        target = re.compile(res_enz['sapi'][rs])
        cut_site = target.search(sequence)
        if cut_site is not None:
            sapi = 'Error: Present'
            break
        else: 
            sapi = None
    seq_results.append(sapi)

    # Test for NotI restriction site. Must be absent
    target = re.compile(res_enz['noti'])
    cut_site = target.search(sequence)
    if cut_site is not None:
        noti = 'Error: Present'
    else: 
        noti = None
    seq_results.append(boti)


    # Test for NdeI restrictions site. Must be presenet once
    # upstream and once downstream
    target = re.compile(res_enz['ndei']['ndei_regex'])
    if target.search(sequence) is not None:
        ndei = None

        seq_cut1 = sequence.count(res_enz['ndei']['ndei_cut_f'])
        seq_cut2 = sequence.count(res_enz['ndei']['ndei_cut_r'])
        if seq_cut1 != 1:
            ndei = 'Error in NdeI cut site'
        if seq_cut2 != 1:
            ndei = 'Error in NdeI cut site'
    seq_results.append(ndei)

    if bsai and sapi and noti and ndei is not None:
        return seq_results
    else:
        return None

    return seq_results

def get_primer_pairs(f_primers, r_primers):
    name = []
    sequence = []

    with open(f_primers) as forward, open(r_primers) as reverse:
        for line1, line2 in zip(forward, reverse):
            if '>' in line1:
                name.append((line1.rsplit(), line2.rsplit()))
            else:
                sequence.append((str(line1.rsplit()), str(line2.rsplit())))
    print sequence

            

def primer_test(primer_pairs, sequence):

    # List to hold sequence of unused pairs
    unused_pairs = []
    
    for pair in primer_pairs:
        # Regex to search for sequence which
        # can be amplified using the primer-pair    
        regex = '^.*{forward}.+{reverse}.*$'.format(
            forward=pair[0], reverse=[pair[1]])
        target = re.compile(regex)
        
        # Uses grep to find matches in the sequence list.
        # Returns the number of matches
        count = subprocess.Popen(['grep', '-c', regex], stdout=subprocess.PIPE)
        working_pair = int(count.communicate())
        
        # Appends an unused primer pair to list of unused primers.
        if working_pair == 0:
            unused pair.append([pair])


if __name__ == "__main__":

    # Sets filetype for t)he what is being tested
    filetype = 'Primer'
    filetype = 'Promoter'

    primer_pairs = get_primer_pairs(f_primers, r_primers)

    with open(inputfile) as fh:
        for fasta in SeqIO.parse(inputfile, 'fasta'):
            name = fasta.id
            sequence = fasta.seq.upper()
            
            # Test the length of every sequence
            length_errors.append(length_test(name, sequence, filetype))

            if filetype == 'Promoter':
                # Test for restriction sites
                res_results = res_site_test(name, sequence, res_enz)
                if res_results is not None:
                    res_errors.append(res_results)

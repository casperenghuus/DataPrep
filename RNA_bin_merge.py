import os
import re

from tabulate import tabulate

import NGS_CNE
import RNA_contigs

CWD = '/scratch/cne/ecre/counts/202_hsrna/'
# CWD = os.getcwd()

if __name__ == "__main__":

    # Sets input-type (RNA or DNA)
    ftype = 'RNA' #CHANGE!
    if ftype == 'DNA':
        all_results = os.path.join(CWD, 'sequences.allresults.dna.txt')
        all_stats = os.path.join(CWD, 'stats.allresults.dna.txt')
    if ftype == 'RNA':
        all_results = os.path.join(CWD, 'sequences.allresults.rna.txt')
        all_stats = os.path.join(CWD, 'stats.allresults.rna.txt')
    if ftype == 'Small':
        all_results = os.path.join(CWD, 'sequences.allresults.small.rna.txt')
        all_stats = os.path.join(CWD, 'stats.allresults.small.rna.txt')

    if ftype == 'Small':
        RNA_small_contigs.seq_counts()
        RNA_small_contigs.merge_bins(ftype)
        RNA_small_contigs.run_bowtie()
        RNA_small_contigs.bowtie_to_file(all_results)
        RNA_small_contigs.final_stats(all_stats, all_results)

    if ftype != 'Small':
        # RNA_contigs.seq_counts()
        # RNA_contigs.merge_bins(ftype)
        # RNA_contigs.run_bowtie()
        RNA_contigs.bowtie_to_file(all_results)
        RNA_contigs.final_stats(all_stats, all_results)



import subprocess
from NGS_CNE import load_fq_files


if __name__ == '__main__':

	file_pairs = load_fq_files()

	for i in file_pairs:
		print 'python run_seqprep.py {forward_file} {reverse_file}'.format(
			forward_file=i[0], reverse_file=i[1])

# BASH COMMAND:
# python run_all_pairs.py | parallel -P 18 &
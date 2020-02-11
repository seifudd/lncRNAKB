import numpy as np
import pandas as pd 
import os
import argparse

parser = argparse.ArgumentParser(description='Extract important metrics from Alignment data results.')
parser.add_argument('--tissue_name', type=str)

def alignment_data(args):
	current_dir = os.getcwd() 
	data_directory = os.path.join(current_dir, args.tissue_name)
	csv_folder = os.path.join(current_dir, 'CSV_Files')
	data_dict = {}
	for filename in os.listdir(data_directory):
		if filename.endswith('.slurm.err.txt'):
			run_s = filename.split('.')[0]
			data_dict[run_s] = []
			with open(os.path.join(data_directory, filename), 'r') as file:
				contents = file.readlines()
				for line in contents:
					if 'started analysis' in line.lower() and '1P.fastq.gz' in line:
						data_dict[run_s].append(line.split()[3])
					if 'started analysis' in line.lower() and '2P.fastq.gz' in line:
						data_dict[run_s].append(line.split()[3])
					if 'reads; of these' in line.lower():
						data_dict[run_s].append(line.split()[0])
					if 'overall alignment rate' in line.lower():
						percent = line.split()[0]
						#percent = ''.join(c for c in percent if c not in '()')
						data_dict[run_s].append(percent)
					if 'concordantly exactly 1 time' in line.lower():
						this_line = line.split()
						num_reads = this_line[0]
						data_dict[run_s].append(num_reads)
						percent = this_line[1]
						percent = ''.join(c for c in percent if c not in '()')
						data_dict[run_s].append(percent)

	tissue_df = pd.DataFrame.from_dict(data_dict, orient='index', columns=['Read1File', 'Read2File', 'PairedEndReadNum', 'UniquePairedEndReadNum', 'UniqueAlignmentRate', 'OverallAlignmentRate'])
	tissue_df.index.name = 'Run_s'

	tissue_df.to_csv(os.path.join(csv_folder, args.tissue_name+'_stats_by_sample.csv'))


if __name__ == '__main__':
	args = parser.parse_args()
	alignment_data(args)
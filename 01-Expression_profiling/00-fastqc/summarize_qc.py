import numpy as np
import pandas as pd 
import os
import argparse
import zipfile
from io import StringIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Summarize FastQC Data')
parser.add_argument('--tissue_name', type=str)


def summarize_qc(args):
	current_dir = os.getcwd() 
	data_directory = os.path.join(current_dir, args.tissue_name)
	csv_folder = os.path.join(current_dir, 'CSV_Files')
	plots_folder = os.path.join(current_dir, 'Plots')
	data = pd.DataFrame()

	for sample_file in os.listdir(data_directory):
		if not sample_file.startswith('.'):
			# Code to unzip the files.
			for filename in os.listdir(os.path.join(data_directory, sample_file)):
				
				if filename.endswith('fastqc.zip'):
					with zipfile.ZipFile(os.path.join(data_directory,sample_file,filename), 'r') as zip_ref:
						unzipped_name = filename.split('.')[0]
						zip_ref.extractall(os.path.join(data_directory,sample_file))


			# Gets quality data from Read 1.
			with open(os.path.join(data_directory, sample_file, sample_file+'_1P_fastqc','fastqc_data.txt')) as file:
				content = file.read().replace('#', '>>').split('>>')
		
				sequence_quality_table = StringIO(content[7])
				sequence_quality_df_read1 = pd.read_table(sequence_quality_table)


			# Gets quality data from Read 2.
			with open(os.path.join(data_directory, sample_file, sample_file+'_2P_fastqc','fastqc_data.txt')) as file:
				content = file.read().replace('#', '>>').split('>>')
		
				sequence_quality_table = StringIO(content[7])
				sequence_quality_df_read2 = pd.read_table(sequence_quality_table)


			# Take means for each read and concatenate.
			sequence_quality_mean_concat = pd.concat([sequence_quality_df_read1[['Mean']], sequence_quality_df_read2[['Mean']]], axis=1)
			sequence_quality_mean_concat.columns = ['read1_'+sample_file, 'read2_'+sample_file]

			if data.empty:
				data = sequence_quality_df_read1[['Base']]

			# Put into one file.
			data = pd.concat((data, sequence_quality_mean_concat), axis=1)


	read1_data = data.filter(like='read1', axis=1)
	read2_data = data.filter(like='read2', axis=1)


	plt.figure()
	ax = read1_data.plot(legend=False)
	plt.suptitle('Read 1 Quality Scores')
	plt.ylabel('Quality')
	plt.xlabel('Base')
	plt.savefig(os.path.join(plots_folder, args.tissue_name+'read1_plot.png'))

	plt.figure()
	ax = read2_data.plot(legend=False)
	plt.suptitle('Read2 Quality Scores')
	plt.xlabel('Base')
	plt.ylabel('Quality')
	plt.savefig(os.path.join(plots_folder, args.tissue_name+'read2_plot.png'))

	read1_data = pd.concat((data[['Base']], read1_data), axis=1)
	read2_data = pd.concat((data[['Base']], read2_data), axis=1)

	data.to_csv(os.path.join(csv_folder, args.tissue_name+'quality_spreadsheet.csv'))

	# read1_data.to_csv(os.path.join(csv_folder, args.tissue_name+'read1_quality_spreadsheet.csv'))
	# read2_data.to_csv(os.path.join(csv_folder, args.tissue_name+'read2_quality_spreadsheet.csv'))


# plt.figure()
# data.T.plot(legend=False)
# plt.savefig('plot.png')

if __name__ == '__main__':
	args = parser.parse_args()
	summarize_qc(args)

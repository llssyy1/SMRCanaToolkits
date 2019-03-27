

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import copy
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
plt.switch_backend('agg');

description = "";
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=True);
parser.add_argument("-n", "--names", required=True, help="names of samples");
##
# input
# SAMPLES
# output:
# sample+"_Events/", barh, bar,	barh_align
# args = parser.parse_args();
# SAMPLES= args.names.replace(' ', '').split(",");
#SAMPLES=["K140", "K510", "SEC", "SHEE", "TE5"];

###############

def main():
	args = parser.parse_args();
	SAMPLES= args.names.replace(' ', '').split(",");
	# run suppa generateEvents to get ioe and ioi file
	for sample in SAMPLES:
		gtf=sample+".gtf";
		directory=sample+"_Events/";
		if(not os.path.exists(directory)):
			os.mkdir(directory);
		header=directory+sample;
		os.system("python ~/miniconda2/envs/test_py3/bin/suppa.py generateEvents -i "+gtf+" -o "+header+" -e SE SS MX RI FL -f ioe");
		os.system("python ~/miniconda2/envs/test_py3/bin/suppa.py generateEvents -i "+gtf+" -o "+header+" -f ioi");

	def plot_bar(df):
		n_row=df.shape[0];
		n_col=df.shape[1];
		df.plot(kind='bar', figsize=(2*n_col,12));# 横柱形图 kind=barh
		plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., markerscale=10);
		plt.savefig('bar.png', dpi=600, format='png');
		plt.show();

	def plot_barh(df):
		n_row=df.shape[0];
		n_col=df.shape[1];
		df.plot(kind='bar', stacked=True, figsize=(2*n_col,12));# stacked=True堆叠
		plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., markerscale=10);
		plt.savefig('barh.png', dpi=600, format='png');
		plt.show();

	def plot_barh_align(df):
		mat_proportion=copy.copy(df);
		sum_mat=df.sum(1).tolist();
		n_row=df.shape[0];
		n_col=df.shape[1];
		for i in range(0, n_row):
			for j in range(0, n_col):
				mat_proportion.iloc[i, j]=mat_proportion.iloc[i, j]/sum_mat[i];
		# mat_proportion.plot(kind='bar', stacked=True, figsize=(20,12));
		mat_proportion.plot(kind='bar', stacked=True, figsize=(2*n_col,12));
		plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., markerscale=10);
		plt.savefig('barh_align.png', dpi=600, format='png');
		plt.show();

	##
	def plot_AS(SAMPLES):
		data_DF=[];
		ASs=["A3", "A5", "AF", "AL", "MX", "RI", "SE"];
		for sample in SAMPLES:
			labels = [];
			fracs = [];
			for AS in ASs:
				directory=sample+"_Events/";
				IOE=directory+sample+"_"+AS+"_strict.ioe";
				count = len(open(IOE, 'r').readlines())-1;
				labels.append(AS);
				fracs.append(count);
			print(sample+":");
			print(labels);
			print(fracs);
			data_DF.append(fracs);
		mat=pd.DataFrame(data=data_DF, index=SAMPLES, columns=ASs);
		mat.to_csv("AS sta.txt", sep="\t");
		plot_barh(mat);
		plot_bar(mat);
		plot_barh_align(mat);

	plot_AS(SAMPLES);


if __name__ == '__main__':
    main()





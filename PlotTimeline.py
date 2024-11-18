#!/usr/bin/env python

# (25/10/2024)	Script to plot a timeline of the evolutionary dynamics using the output data of ProtEvol.

# Import modules #

import math, sys, os
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=0, gap_score=0)

CodonTable = {'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L','TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C','TGA':'*','TGG':'W','CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P','CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R','ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T','AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R','GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G'}

# Define functions #

def PrintHelp():
	print("Error --- usage: ./PlotTimeline.py -d [data file] -i [input sequence] -T [target structure] -o [figure name]")
	sys.exit(1)

def Translate(nt_seq):
	S = ""
	for i in range(0, len(nt_seq)-len(nt_seq)%3,3):
		S += CodonTable[nt_seq[i:i+3]]
	return S.split('*')[0]

def AlignStrings(string1, string2):
	if len(string1) == 0 or len(string2) == 0:
		return 0
	else:
		return aligner.score(string1, string2)
	
# Parse arguments #

data_file = ""
input_sequence = ""
target_structure = ""
figure_name = ""

i_arg = 1
while i_arg != len(sys.argv):
	if sys.argv[i_arg] == '-d':
		data_file = sys.argv[i_arg+1]
	elif sys.argv[i_arg] == '-i':
		input_sequence = sys.argv[i_arg+1]
	elif sys.argv[i_arg] == '-T':
		target_structure = sys.argv[i_arg+1]
	elif sys.argv[i_arg] == '-o':
		figure_name = sys.argv[i_arg+1]
	else:
		PrintHelp()
	i_arg += 2

if data_file == "" or input_sequence == "" or target_structure == "" or figure_name == "":
	PrintHelp()

# Load data #

Data = pd.read_csv(data_file, sep='\t', names=['time','idx','parent_idx','nt_sequence','aa_sequence','ss_structure','exposure','stability','complexity','fitness'])

# Manipulate data #

for col in ['nt_sequence','aa_sequence','ss_structure','exposure']:
	Data.loc[Data[col].isna(),col] = ""


#Assign sequence length
Data['sequence_length'] = Data['aa_sequence'].str.len()

#Assign similarity to target
Data['target_struct_similarity'] = Data['ss_structure'].apply(lambda x: AlignStrings(x, target_structure))

#Assign similarity to parent
Data['parent_struct_similarity'] = 0.
Data.loc[Data['time']==0, 'parent_struct_similarity'] = Data[Data['time']==0]['ss_structure'].str.len()
Data.loc[Data['time']!=0, 'parent_struct_similarity'] = Data[Data['time']!=0].apply(lambda x: AlignStrings(Data[(Data['idx']==x['parent_idx']) & (Data['time']==x['time']-1)]['ss_structure'].values[0], x['ss_structure']), axis=1)

#Assign sequence similarity to input sequence
Data['initial_seq_similarity'] = Data['aa_sequence'].apply(lambda x: AlignStrings(x, input_sequence))

MeanData = Data.groupby(by='time').agg('mean', numeric_only=True).reset_index()

# Plot #

sns.set_theme()
sns.set_style("ticks", {"axes.grid": True, "grid.color": "lightgrey", "grid.linestyle": ":"})
sns.set_context("notebook", font_scale=0.7)

fig, axs = plt.subplots(6, 1, figsize=(10,12), tight_layout=True, sharex=True)

plot_vars = ['sequence_length','initial_seq_similarity','parent_struct_similarity','target_struct_similarity', 'fitness', 'complexity']
plot_colors = ['black', 'brown', 'gold', 'midnightblue', 'green', 'pink']

axs[0].set_xlim([0,1000])

for i, var in enumerate(plot_vars):
	if var not in MeanData.columns:	continue
	axs[i].plot(Data['time'], Data[var], 's', color=plot_colors[i], ms=6, clip_on=False, label="", alpha=0.03, mew=0)
	axs[i].plot(MeanData['time'], MeanData[var], '-', color=plot_colors[i])
	if var not in ['fitness','complexity']:	axs[i].set_ylim([0,300])
	axs[i].set_ylabel(var)
	axs[i].grid(True)

plt.savefig(figure_name)
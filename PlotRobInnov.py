#!/usr/bin/env python

# (15/06/2026)	Visualize robustness to mutations and innovation across various mutation rates for different proteins.

# Import modules #

import math, sys, os
from Bio.Align import PairwiseAligner
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=0, gap_score=0)

def ReadSS(filename):
	Dict = {'Wildtype':[],'Sequence':[],'Secondary structure':[]}
	cntr = 0
	with open(filename, 'r') as fin:
		for line in fin:
			line = line.rstrip()
			if line[0] == '>':
				lab = line[1:].split(' ')[0]
				cntr = 0
			elif cntr == 1:
				seq = line
			elif cntr == 2:
				Dict['Wildtype'].append(lab)
				Dict['Sequence'].append(seq)
				Dict['Secondary structure'].append(line.split(' ')[0])
			cntr += 1
	return Dict

def AlignStrings(string1, string2):
	if len(string1) == 0 or len(string2) == 0:
		return 0
	else:
		return aligner.score(string1, string2)

####################
### Protein data ###
####################

computers = ["binfgpu8", "binfgpu8", "binfgpu5", "binfgpu5"]
projects = ["III","IV","V","VI"]
mut_rates = [0.01, 0.1, 0.001, 0.0001]
prot_wildtypes = ["A0A178VEK7", "A5LHX3", "O14842", "O15144", "O35980", "O88792", "O94457", "O95407", "O95749", "P08191"]
columns = ['time','idx','parent_idx','nt_sequence','aa_sequence','ss_structure','exposure','stability','complexity','fitness']

#Gather data.
frames = []
for p, m, c in zip(projects, mut_rates, computers):
	for i in range(1,11):
		D = pd.read_csv(f"/hosts/linuxhome/{c}/tmp/sam/protevol/Nimwegen{p}{i}/Nimwegen{p}{i}.out", sep="\t", names=columns)
		D["mutation_rate"] = m
		D["wildtype"] = prot_wildtypes[i-1]
		D = D[D["time"].isin([0,1000])]
		frames.append(D)
		
Data = pd.concat(frames, ignore_index=True, axis=0)
Data['initial_seq'] = Data.apply(lambda x: Data[(Data["time"]==0) & (Data["wildtype"]==x["wildtype"])]["aa_sequence"].values[0], axis=1)
Data['initial_seq_similarity'] = Data.apply(lambda x: AlignStrings(x["aa_sequence"], x["initial_seq"]), axis=1)
Data = Data[Data["time"]==1000]

Data["Structure divergence"] = 300 - Data["fitness"]
Data["Sequence divergence"] = 300 - Data["initial_seq_similarity"]

DM = Data.groupby(by=["time","mutation_rate","wildtype"], as_index=False).min()
DM = DM.sort_values(by="mutation_rate")
print(DM)


################
### RNA data ###
################

rna_wildtypes = ["URS00000034F0", "URS0000005C90", "URS00000114AA", "URS00000200E9", "URS0000025426", "URS0000042CD1", "URS000004AEB1", "URS000004C426", "URS000005E030", "URS00000644DC"]
columns = ['time','idx','parent_idx','sequence','structure','fitness','complexity']
mut_rates = [0.0001, 0.001, 0.01, 0.1]

#Gather data.
frames = []
for i, m in enumerate(mut_rates):
	for k in range(1,11):
		R = pd.read_csv(f"/net/vacuole1/linuxhome/ph-group/sam/NOBINFBACKUP-Documents/ProteinEvol/Yeast/MyAnalysis/alternative_folding/simulation/rna_robu_inno/robu_inno2_m{i}_i{k}.dat", sep="\t", names=columns)
		R["mutation_rate"] = m
		R["wildtype"] = rna_wildtypes[k-1]
		R = R[R["time"].isin([0,1000])]
		frames.append(R)

#Read initial seqs from separate file.
RD = ReadSS("/net/vacuole1/linuxhome/ph-group/sam/NOBINFBACKUP-Documents/ProteinEvol/Yeast/MyAnalysis/alternative_folding/simulation/rna_robu_inno/L300_rnacentral_10wts.ss")

Rata = pd.concat(frames, ignore_index=True, axis=0)
Rata['initial_seq'] = Rata['wildtype'].apply(lambda x: RD['Sequence'][RD['Wildtype'].index(x)])
Rata['initial_seq_similarity'] = Rata.apply(lambda x: AlignStrings(x["sequence"], x["initial_seq"]), axis=1)
Rata = Rata[Rata["time"]==1000]

Rata["Structure divergence"] = 300 - Rata["fitness"]
Rata["Sequence divergence"] = 300 - Rata["initial_seq_similarity"]

RM = Rata.groupby(by=["time","mutation_rate","wildtype"], as_index=False).min()
RM = RM.sort_values(by="mutation_rate")
print(RM)

################
### Plotting ###
################

sns.set_theme()
sns.set_style("ticks", {"axes.grid": True, "grid.color": "lightgrey", "grid.linestyle": ":"})
sns.set_context("talk")

fig, axs = plt.subplots(2, 3, figsize=(15,10), tight_layout=True)



sns.scatterplot(ax=axs[0][0], data=DM[DM["time"]==1000], x="mutation_rate", y="Structure divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
for wt in prot_wildtypes:
	axs[0][0].plot(DM[DM["wildtype"]==wt]["mutation_rate"].values, DM[DM["wildtype"]==wt]["Structure divergence"].values, 'k-', alpha=0.3, zorder=2)

axs[0][0].set_xscale("log")
axs[0][0].set_ylim([0,300])


sns.scatterplot(ax=axs[0][1], data=DM[DM["time"]==1000], x="mutation_rate", y="Sequence divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
for wt in prot_wildtypes:
	axs[0][1].plot(DM[DM["wildtype"]==wt]["mutation_rate"].values, DM[DM["wildtype"]==wt]["Sequence divergence"].values, 'k-', alpha=0.3, zorder=2)

axs[0][1].set_xscale("log")
axs[0][1].set_ylim([0,300])


sns.scatterplot(ax=axs[0][2], data=DM[DM["time"]==1000], x="Structure divergence", y="Sequence divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
for wt in prot_wildtypes:
	axs[0][2].plot(DM[DM["wildtype"]==wt]["Structure divergence"].values, DM[DM["wildtype"]==wt]["Sequence divergence"].values, 'k-', alpha=0.3, zorder=2)

axs[0][2].set_xlim([0,300])
axs[0][2].set_ylim([0,300])



sns.scatterplot(ax=axs[1][0], data=RM[RM["time"]==1000], x="mutation_rate", y="Structure divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
for wt in rna_wildtypes:
	axs[1][0].plot(RM[RM["wildtype"]==wt]["mutation_rate"].values, RM[RM["wildtype"]==wt]["Structure divergence"].values, 'k-', alpha=0.3, zorder=2)

axs[1][0].set_xscale("log")
axs[1][0].set_ylim([0,300])


sns.scatterplot(ax=axs[1][1], data=RM[RM["time"]==1000], x="mutation_rate", y="Sequence divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
for wt in rna_wildtypes:
	axs[1][1].plot(RM[RM["wildtype"]==wt]["mutation_rate"].values, RM[RM["wildtype"]==wt]["Sequence divergence"].values, 'k-', alpha=0.3, zorder=2)

axs[1][1].set_xscale("log")
axs[1][1].set_ylim([0,300])


sns.scatterplot(ax=axs[1][2], data=RM[RM["time"]==1000], x="Structure divergence", y="Sequence divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
for wt in rna_wildtypes:
	axs[1][2].plot(RM[RM["wildtype"]==wt]["Structure divergence"].values, RM[RM["wildtype"]==wt]["Sequence divergence"].values, 'k-', alpha=0.3, zorder=2)

axs[1][2].set_xlim([0,300])
axs[1][2].set_ylim([0,300])



plt.savefig("RobuInno3.png")
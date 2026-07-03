#!/usr/bin/env python

# (15/06/2026)	Visualize robustness to mutations and innovation across various mutation rates for different proteins.

# Import modules #

import math, sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

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
	sim = 0
	for p in range(len(string1)):
		if string1[p] == string2[p]:
			sim += 1
	return sim

####################
### Protein data ###
####################

computers = ["binfgpu5", "binfgpu5", "binfgpu8", "binfgpu8", "binfgpu2", "binfgpu2"]
projects = ["I","II","III","IV","V","VI"]
mut_rates = [0.001, 0.002, 0.003, 0.005, 0.007, 0.010]
prot_wildtypes = ["A0A178VEK7", "A5LHX3", "O14842", "O15144", "O35980", "O88792", "O94457", "O95407", "O95749", "P08191"]
columns = ['time','idx','parent_idx','nt_sequence','aa_sequence','ss_structure','exposure','stability','complexity','fitness']

#Gather data.
struct_frames = []
frames = []
for p, m, c in zip(projects, mut_rates, computers):
	for i in range(1,11):
		if not os.path.exists(f"/hosts/linuxhome/{c}/tmp/sam/protevol/Huynen{p}{i}/Huynen{p}{i}.out"):
			continue
		D = pd.read_csv(f"/hosts/linuxhome/{c}/tmp/sam/protevol/Huynen{p}{i}/Huynen{p}{i}.out", sep="\t", names=columns)
		if 1000 not in D['time'].values:
			continue
		D["mutation_rate"] = m
		D["wildtype"] = prot_wildtypes[i-1]
		target_structure = D[D['time'] == 0]['ss_structure'].values[0]
		initial_sequence = D[D['time'] == 0]['aa_sequence'].values[0]

		# --- novel structure discovery ---
		S = (D[["wildtype", "time", "aa_sequence", "ss_structure"]].loc[D['ss_structure'] != target_structure].drop_duplicates(subset='ss_structure', keep='first').sort_values('time').reset_index(drop=True))
		S['ss_distance'] = 300 - S['ss_structure'].apply(lambda x: AlignStrings(x, target_structure))
		S['seq_distance'] = 300 - S['aa_sequence'].apply(lambda x: AlignStrings(x, initial_sequence))
		S['cumulative_unique'] = S.index + 1		  # already sorted by time, one row per novel structure
		S['running_max_dist'] = S['ss_distance'].cummax()
		S['running_max_seqdist'] = S['seq_distance'].cummax()
		S['mutation_rate'] = m
		struct_frames.append(S)

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
# print(DM)

Sprot = pd.concat(struct_frames, ignore_index=True, axis=0)
# print(Sprot)

################
### RNA data ###
################

rna_wildtypes = ["URS00000034F0", "URS0000005C90", "URS00000114AA", "URS00000200E9", "URS0000025426", "URS0000042CD1", "URS000004AEB1", "URS000004C426", "URS000005E030", "URS00000644DC"]
columns = ['time','idx','parent_idx','sequence','structure','fitness','complexity']
# mut_rates = [0.0001, 0.001, 0.01, 0.1]
mut_rates = [0.001, 0.002, 0.003, 0.005, 0.007, 0.01]

#Gather data.
struct_frames = []
frames = []
for i, m in enumerate(mut_rates):
	for k in range(1,11):
		R = pd.read_csv(f"/net/vacuole1/linuxhome/ph-group/sam/NOBINFBACKUP-Documents/ProteinEvol/Yeast/MyAnalysis/alternative_folding/simulation/rna_robu_inno/robu_inno4_m{i}_i{k}.dat", sep="\t", names=columns)
		R["mutation_rate"] = m
		R["wildtype"] = rna_wildtypes[k-1]

		target_structure = R[R['time'] == 0]['structure'].values[0]
		initial_sequence = R[R['time'] == 0]['sequence'].values[0]

		# --- novel structure discovery ---
		S = (R[["wildtype", "time", "sequence", "structure"]].loc[R['structure'] != target_structure].drop_duplicates(subset='structure', keep='first').sort_values('time').reset_index(drop=True))
		S['struct_distance'] = 300 - S['structure'].apply(lambda x: AlignStrings(x, target_structure))
		S['seq_distance'] = 300 - S['sequence'].apply(lambda x: AlignStrings(x, initial_sequence))
		S['cumulative_unique'] = S.index + 1		  # already sorted by time, one row per novel structure
		S['running_max_dist'] = S['struct_distance'].cummax()
		S['running_max_seqdist'] = S['seq_distance'].cummax()
		S['mutation_rate'] = m
		struct_frames.append(S)

		R = R[R["time"].isin([0,1000])]
		frames.append(R)

Rata = pd.concat(frames, ignore_index=True, axis=0)
Rata['initial_seq'] = Rata.apply(lambda x: Rata[(Rata["time"]==0) & (Rata["wildtype"]==x["wildtype"])]["sequence"].values[0], axis=1)
Rata['initial_seq_similarity'] = Rata.apply(lambda x: AlignStrings(x["sequence"], x["initial_seq"]), axis=1)
Rata = Rata[Rata["time"]==1000]

Rata["Structure divergence"] = 300 - Rata["fitness"]
Rata["Sequence divergence"] = 300 - Rata["initial_seq_similarity"]

RM = Rata.groupby(by=["time","mutation_rate","wildtype"], as_index=False).min()
RM = RM.sort_values(by="mutation_rate")
print(RM)

Srna = pd.concat(struct_frames, ignore_index=True, axis=0)
print(Srna)

################
### Plotting ###
################

# sns.set_theme()
# sns.set_style("ticks", {"axes.grid": True, "grid.color": "lightgrey", "grid.linestyle": ":"})
# sns.set_context("talk")

# fig, axs = plt.subplots(2, 3, figsize=(15,10), tight_layout=True)



# sns.scatterplot(ax=axs[0][0], data=DM[DM["time"]==1000], x="mutation_rate", y="Structure divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
# for wt in prot_wildtypes:
# 	axs[0][0].plot(DM[DM["wildtype"]==wt]["mutation_rate"].values, DM[DM["wildtype"]==wt]["Structure divergence"].values, '-', alpha=1.0, zorder=2)

# axs[0][0].set_xscale("log")
# axs[0][0].set_ylim([0,300])

# DM[r'$\mu\times L\times N\times t$'] = DM['mutation_rate']*10*DM['time']*300

# sns.scatterplot(ax=axs[0][1], data=DM[DM["time"]==1000], x=r'$\mu\times L\times N\times t$', y="Sequence divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
# for wt in prot_wildtypes:
# 	axs[0][1].plot(DM[DM["wildtype"]==wt][r'$\mu\times L\times N\times t$'].values, DM[DM["wildtype"]==wt]["Sequence divergence"].values, '-', alpha=1.0, zorder=2)

# # axs[0][1].set_xscale("log")
# axs[0][1].set_xlim([0,30000])
# axs[0][1].set_ylim([0,300])

# for (wt, m), grp in Sprot.groupby(['wildtype', 'mutation_rate']):
# 	if m!=0.002: continue
# 	# axs[0][2].plot(grp['time'], grp['cumulative_unique'])
# 	axs[0][2].plot(grp['running_max_seqdist'], grp['running_max_dist'])

# axs[0][2].set_xlim([0,300])
# axs[0][2].set_ylim([0,300])
# axs[0][2].set_xlabel("Max. sequence distance")
# axs[0][2].set_ylabel("Max. structure distance")



# sns.scatterplot(ax=axs[1][0], data=RM[RM["time"]==1000], x="mutation_rate", y="Structure divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
# for wt in rna_wildtypes:
# 	axs[1][0].plot(RM[RM["wildtype"]==wt]["mutation_rate"].values, RM[RM["wildtype"]==wt]["Structure divergence"].values, '-', alpha=0.3, zorder=2)

# axs[1][0].set_xscale("log")
# axs[1][0].set_ylim([0,300])

# RM[r'$\mu\times L\times N\times t$'] = RM['mutation_rate']*10*RM['time']*300

# sns.scatterplot(ax=axs[1][1], data=RM[RM["time"]==1000], x=r'$\mu\times L\times N\times t$', y="Sequence divergence", marker='o', lw=0., alpha=1, s=50, clip_on=False, legend=False, color="#66526E", zorder=3)
# for wt in rna_wildtypes:
# 	axs[1][1].plot(RM[RM["wildtype"]==wt][r'$\mu\times L\times N\times t$'].values, RM[RM["wildtype"]==wt]["Sequence divergence"].values, '-', alpha=0.3, zorder=2)

# # axs[1][1].set_xscale("log")
# axs[0][1].set_xlim([0,30000])
# axs[1][1].set_ylim([0,300])

# for (wt, m), grp in Srna.groupby(['wildtype', 'mutation_rate']):
# 	if m!=0.002: continue
# 	# axs[0][2].plot(grp['time'], grp['cumulative_unique'])
# 	axs[1][2].plot(grp['running_max_seqdist'], grp['running_max_dist'])

# axs[1][2].set_xlim([0,300])
# axs[1][2].set_ylim([0,300])
# axs[1][2].set_xlabel("Max. sequence distance")
# axs[1][2].set_ylabel("Max. structure distance")

# plt.savefig("RobuInno8.png")





##################
### Final Plot ###
###################

sns.set_theme()
sns.set_style("ticks", {"axes.grid": True, "grid.color": "lightgrey", "grid.linestyle": ":"})
sns.set_context("talk")

fig, axs = plt.subplots(1, 2, figsize=(10,5), tight_layout=True)

for wt in rna_wildtypes:
	axs[0].plot([x for x in RM[RM["wildtype"]==wt]["mutation_rate"].values], RM[RM["wildtype"]==wt]["Structure divergence"].values, 'o-', lw=1.5, ms=5, alpha=1.0, zorder=2, color='#D40055')

for wt in prot_wildtypes:
	axs[0].plot([x for x in DM[DM["wildtype"]==wt]["mutation_rate"].values], DM[DM["wildtype"]==wt]["Structure divergence"].values, 'o-', lw=1.5, ms=5, alpha=1.0, zorder=2, color='#00AAD4')

axs[0].set_xscale("log")
# axs[0].set_xlim([0.001,0.01])
axs[0].set_ylim([0,150])
axs[0].set_xlabel('Mutation rate')
axs[0].set_ylabel('Distance to target structure')

for (wt, m), grp in Srna.groupby(['wildtype', 'mutation_rate']):
	if m!=0.002: continue
	axs[1].plot(grp['running_max_seqdist'], grp['running_max_dist'], '-', color='#D40055', lw=1.5)

for (wt, m), grp in Sprot.groupby(['wildtype', 'mutation_rate']):
	if m!=0.002: continue
	axs[1].plot(grp['running_max_seqdist'], grp['running_max_dist'], '-', color='#00AAD4', lw=1.5)

axs[1].set_xlim([0,300])
axs[1].set_ylim([0,300])
axs[1].set_xlabel("Most distant sequence")
axs[1].set_ylabel("Most distant structure")


plt.savefig("RobuInnoFinal2.svg")
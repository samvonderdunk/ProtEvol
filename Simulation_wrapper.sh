#!/usr/bin/env bash

# (14/10/2024)	Wrapper for ProtEvol.

# #Thioredoxin

# GPU=1
# ProtEvol=./Simulation.py
# ProjectBaseName=ThioredoxinD
# PopSize=100
# SimTime=1000
# CompScale=2
# StartSeed=315

# SimDir=/net/vacuole1/linuxhome/ph-group/sam/NOBINFBACKUP-Documents/ProteinEvol/Yeast/MyAnalysis/alternative_folding/simulation/proteins
# TargetStruct=`cat ${SimDir}/Q8NBQ5.sss`
# InitSeq=`cat ${SimDir}/Q55C17.seq`

# i=1
# for SelectTargetP in 0.1 0.25 0.5 0.75 0.9;
# do
# 	for MutRate in 0.0001 0.00025 0.0005;
# 	do

# 		Seed=$((StartSeed+i))
# 		ProjectName=${ProjectBaseName}${i}

# 		CUDA_VISIBLE_DEVICES=$GPU $ProtEvol -p $ProjectName -s $Seed -m $MutRate -N $PopSize -t $SimTime -C $CompScale -T $TargetStruct -i $InitSeq -pst $SelectTargetP > /linuxhome/tmp/sam/protevol/${ProjectName}.log

# 		i=$((i+1))

# 	done
# done

#RPS
ProtEvol=/mnt/users/dunks/Projects/ProteinEvolution/ProtEvol/Simulation.py
ProjectBaseName=RPS
PopSize=100
SimTime=1000
CompScale=0
StartSeed=330
MutRate=0.005

SimDir=/mnt/users/dunks/Projects/ProteinEvolution/Simulations/proteins
InitSeqs=(Seq0123.seq Seq0456.seq)

# for i in {41..42};
for i in 41;
do
	Seed=$((StartSeed+i))
	ProjectName=${ProjectBaseName}${i}

	if [ "$i" -lt "36" ]; then
		InitSeq=${InitSeqs[0]}
	else
		InitSeq=${InitSeqs[1]}
	fi

	addqueue -q heraclesgpu -s --gpus 1 --gputype rtx6000adawith48gb -m 25 /usr/local/shared/python/3.9.6/bin/python3 $ProtEvol -p $ProjectName -s $Seed -m $MutRate -N $PopSize -t $SimTime -C $CompScale -i `cat proteins/$InitSeq` -F simplicity_structure -g aa -mu_insertion 0 -mu_deletion 0 -mu_duplication 0 -mu_ablation 0 -mu_reversion 0 -mu_transposition 0 > /mnt/users/dunks/Projects/ProteinEvolution/Simulations/${ProjectName}.log
done
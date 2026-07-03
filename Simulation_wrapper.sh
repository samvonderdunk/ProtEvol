#!/usr/bin/env bash

# (14/10/2024)	Wrapper for ProtEvol.

#Thioredoxin

GPU=1
ProtEvol=./Simulation.py
ProjectBaseName=HuynenVI
PopSize=10
SimTime=1000
CompScale=0
StartSeed=471

SimDir=/net/vacuole1/linuxhome/ph-group/sam/NOBINFBACKUP-Documents/ProteinEvol/Yeast/MyAnalysis/alternative_folding/simulation/proteins
# TargetStruct=`cat ${SimDir}/Q8NBQ5.sss`
# InitSeq=`cat ${SimDir}/Q55C17.seq`
WTS=(A0A178VEK7 A5LHX3 O14842 O15144 O35980 O88792 O94457 O95407 O95749 P08191)

for i in {1..10};
do
	Seed=$((StartSeed+i))
	ProjectName=${ProjectBaseName}${i}
	TargetStruct=`cat ${SimDir}/${WTS[i-1]}.ss3`
	InitSeq=`cat ${SimDir}/${WTS[i-1]}.seq`

	CUDA_VISIBLE_DEVICES=$GPU $ProtEvol -p $ProjectName -s $Seed -m 0.01 -N $PopSize -t $SimTime -C $CompScale -T $TargetStruct -i $InitSeq -mu_insertion 0.0 -mu_deletion 0.0 -mu_duplication 0.0 -mu_ablation 0.0 -mu_reversion 0.0 -mu_transposition 0.0 > /linuxhome/tmp/sam/protevol/${ProjectName}.log
done

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

# #RPS
# GPU=1
# ProtEvol=./Simulation.py
# ProjectBaseName=RPS
# PopSize=100
# SimTime=1000
# CompScale=0
# StartSeed=290
# MutRate=0.005

# SimDir=/net/vacuole1/linuxhome/ph-group/sam/NOBINFBACKUP-Documents/ProteinEvol/Yeast/MyAnalysis/alternative_folding/simulation/proteins
# InitSeq=`cat ${SimDir}/Q55C17.seq`

# for i in {11..20};
# do
# 	Seed=$((StartSeed+i))
# 	ProjectName=${ProjectBaseName}${i}

# 	CUDA_VISIBLE_DEVICES=$GPU $ProtEvol -p $ProjectName -s $Seed -m $MutRate -N $PopSize -t $SimTime -C $CompScale -i $InitSeq -F simplicity_structure -g aa > /linuxhome/tmp/sam/protevol/${ProjectName}.log
# done
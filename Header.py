#!/usr/bin/env python

# Here I just import all modules and define some fixed variables.

import random as rn
import numpy as np
import torch
import math, sys, os
import KC

model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")

import esm
os.environ['MKL_THREADING_LAYER'] = 'GNU'	#Solves an issue with numpy.

model = esm.pretrained.esmfold_v1()
model = model.eval().cuda()
# model.set_chunk_size(128)

linuxhome_dir = '/linuxhome/tmp/sam/protevol'
Nucleotides = ['A','C','G','U']
MutProbsNT = [[0, 1/6, 4/6, 1/6], [1/6, 0, 1/6, 4/6], [4/6, 1/6, 0, 1/6], [1/6, 4/6, 1/6, 0]]
AminoAcids = ['A','G','I','L','P','V','F','W','Y','D','E','R','H','K','S','T','C','M','N','Q']
MutProbsAA = np.ones((20,20))

CodonTable = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

#Table from wikipedia but source is Tien et al. 2013 (empirical data).
MaxExposureTable = {
	'A': 121,
	'R': 265,
	'N': 187,
	'D': 187,
	'C': 140,
	'E': 214,
	'Q': 214,
	'G': 97,
	'H': 216,
	'I': 195,
	'L': 191,
	'K': 230,
	'M': 203,
	'F': 228,
	'P': 154,
	'S': 143,
	'T': 163,
	'W': 264,
	'Y': 255,
	'V': 165
}
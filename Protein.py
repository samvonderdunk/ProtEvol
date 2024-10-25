#!/usr/bin/env python

# Import modules #

from Header import *
import Config

#####################
### Protein class ###
#####################

class Protein:
    
	def __init__(self, idx):
		self.idx = idx
		self.parent_idx = None
		self.nt_sequence = ""
		self.aa_sequence = ""
		self.ss_structure = ""
		self.exposure = ""
		self.stability = None
		self.complexity = None
		self.fitness = None

	def Replicate(self, Parent):
		self.parent_idx = Parent.idx
		self.nt_sequence = Parent.nt_sequence
		self.aa_sequence = Parent.aa_sequence

	def CopyPhenotype(self, Parent):
		self.ss_structure = Parent.ss_structure
		self.exposure = Parent.exposure
		self.stability = Parent.stability
		self.complexity = Parent.complexity
		self.fitness = Parent.fitness

	def Translate(self):
		S = ""
		for i in range(0, len(self.nt_sequence)-len(self.nt_sequence)%3,3):
			S += CodonTable[self.nt_sequence[i:i+3]]
			if S[-1] == '*':
				if rn.random() < Config.transl_stop_prob:
					S = S[:-1]
					break
		self.aa_sequence = S.replace('*','')
		# self.aa_sequence = S.split('*')[0]

	def Mutate(self):

		if Config.genotype_level == 'nt':
			Alphabet = Nucleotides
			MutProbs = MutProbsNT
			s = self.nt_sequence
		elif Config.genotype_level == 'aa':
			Alphabet = AminoAcids
			MutProbs = MutProbsAA
			s = self.aa_sequence

		S = ""
		i = 0
		while i < len(s):
			if rn.random() < Config.mutation_rate:
				mut_type = rn.random()
				if mut_type < Config.p_insertion:	#Insertion.
					S += rn.choice(Alphabet)
					while rn.random() < 0.5:	#Exponentially decaying insertion length.
						S += rn.choice(Alphabet)
					S += s[i]
				elif mut_type < (1.0 - Config.p_deletion):	#Substitution.
					q = rn.choices(Alphabet, weights=MutProbs[Alphabet.index(s[i])], k=1)[0]
					S += q
				else:	#Deletion.
					i += 1
					while rn.random() < 0.5:	#Exponentially decaying deletion length.
						i += 1
			else:
				S += s[i]
				i += 1

		if Config.genotype_level == 'nt':
			self.nt_sequence = S
			self.Translate()

		elif Config.genotype_level == 'aa':
			self.aa_sequence = S

		return S == s	#Report back whether a mutation occurred.

	def ReadSSP(self, ssp_file):
		with open(ssp_file, 'r') as fin:
			read = False
			Sequence = ""
			Structure = ""
			Exposure = ""
			for line in fin:
				line = line.rstrip()
				if read and line[9] != " ":	#For some reason DSSP sometimes inserts exclamation marks, i.e. extra residues in the output, perhaps because the structure does not work, but we just ignore it.
					Structure += line[16]
					Sequence += line[13]
					rel_exp = int(line[35:38])/MaxExposureTable[Sequence[-1]]
					Exposure += ('E' if rel_exp>Config.exp_threshold else 'B')
				if line[:5] == "  #  ": read = True

			if Config.phenotype_level == "1D":
				for ss in ['I','G']:
					Structure = Structure.replace(ss, 'H')
				for ss in [' ','T','S','B']:
					Structure = Structure.replace(ss, 'C')
			elif Config.phenotype_level == "2D":
				Structure = Structure.replace(' ','A')
		return Sequence, Structure, Exposure
	
	def ReadPDB(self, pdb_file):
		b_facts = []
		with open(pdb_file, 'r') as fin:
			for line in fin:
				line = line.rstrip()
				if line[:4] != "ATOM":	continue
				b_facts.append(float(line[61:66]))
		stab = np.mean(b_facts, axis=0)
		if stab == None:
			return 0.0
		else:
			return stab
		
	def BinaryStructure(self):
		s = self.ss_structure
		s = s.replace('A','000')
		s = s.replace('E','111')
		s = s.replace('B','010')
		s = s.replace('T','001')
		s = s.replace('G','011')
		s = s.replace('H','101')
		s = s.replace('I','110')
		s = s.replace('S','100')
		return s
	
	def CalcComplexity(self):
		self.complexity = KC.calc_KC(self.BinaryStructure())

	def MakePhenotype(self):
		if self.aa_sequence == "":
			self.ss_structure = ""
			self.exposure = ""
			self.stability = None
			self.complexity = None

		else:
			with torch.no_grad(), open(f'{linuxhome_dir}/{Config.project_name}/tmp/{self.idx}.pdb', 'w') as fout:
				pdb = model.infer_pdb(self.aa_sequence)
				fout.write(pdb+'\n')

			os.system(f"/home/sam/miniconda3/envs/biolib/bin/mkdssp -i {linuxhome_dir}/{Config.project_name}/tmp/{self.idx}.pdb -o {linuxhome_dir}/{Config.project_name}/tmp/{self.idx}.ssp")

			#Extract secondary structure and exposure strings
			Seq, Str, Exp = self.ReadSSP(f'{linuxhome_dir}/{Config.project_name}/tmp/{self.idx}.ssp')
			self.ss_structure = Str
			self.exposure = Exp

			#Extract stability
			Stab = self.ReadPDB(f'{linuxhome_dir}/{Config.project_name}/tmp/{self.idx}.pdb')
			self.stability = Stab

			#Calculate structural complexity.
			self.CalcComplexity()

	def SimilarityToPhenotype(self, target):
		# # Currently using Hamming distance but could also implement alignment.
		# f = 0
		# for p, q in zip(self.ss_structure, target):
		# 	if p == q:	f += 1
		# return f
	
		if len(self.ss_structure) == 0:
			return 0
		else:
			return aligner.score(self.ss_structure, target)

	def AssignFitness(self):

		if Config.fitness_criterion == 'neutral':
			self.fitness = 100

		elif Config.fitness_criterion == 'target_structure':
			self.fitness = self.SimilarityToPhenotype(Config.target_structure)

		elif Config.fitness_criterion == 'complexity_structure':
			self.fitness = self.complexity

		elif Config.fitness_criterion == 'simplicity_structure':
			self.fitness = - self.complexity

	def PrintString(self):
		return f"{self.idx}\t{self.parent_idx}\t{self.nt_sequence}\t{self.aa_sequence}\t{self.ss_structure}\t{self.exposure}\t{self.stability}\t{self.complexity}\t{self.fitness}"
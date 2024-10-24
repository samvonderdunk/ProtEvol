#!/usr/bin/env python

# Import modules #

from PE5_header import *

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

	def Translate(self):
		S = ""
		for i in range(0, len(self.nt_sequence)-len(self.nt_sequence)%3,3):
			S += CodonTable[self.nt_sequence[i:i+3]]
		self.aa_sequence = S

	def Mutate(self):

		if genotype_level == 'nt':
			Alphabet = Nucleotides
			MutProbs = MutProbsNT
			s = self.nt_sequence
		elif genotype_level == 'aa':
			Alphabet = AminoAcids
			MutProbs = MutProbsAA
			s = self.aa_sequence

		S = ""
		i = 0
		while i < len(s):
			if rn.random() < mutation_rate:
				mut_type = rn.random()
				if mut_type < p_insertion:	#Insertion.
					S += rn.choice(Alphabet)
					while rn.random() < 0.5:	#Exponentially decaying insertion length.
						S += rn.choice(Alphabet)
					S += s[i]
				elif mut_type < (1.0 - p_deletion):	#Substitution.
					q = rn.choices(Alphabet, weights=MutProbs[Alphabet.index(s[i])], k=1)[0]
					S += q
				else:	#Deletion.
					i += 1
					while rn.random() < 0.5:	#Exponentially decaying deletion length.
						i += 1
			else:
				S += s[i]

		if genotype_level == 'nt':
			self.nt_sequence = S
			self.aa_sequence = self.Translate()

		elif genotype_level == 'aa':
			self.aa_sequence = S

	def ReadSSP(ssp_file):
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
					Exposure += ('E' if rel_exp>exp_threshold else 'B')
				if line[:5] == "  #  ": read = True

			if phenotype_level == "1D":
				for ss in ['I','G']:
					Structure = Structure.replace(ss, 'H')
				for ss in [' ','T','S','B']:
					Structure = Structure.replace(ss, 'C')
			elif phenotype_level == "2D":
				Structure = Structure.replace(' ','A')
		return Sequence, Structure, Exposure
	
	def ReadPDB(pdb_file):
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
		s = self.structure
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
		with torch.no_grad():
			print(self.aa_sequence)
			pdb = model.infer_pdb(self.aa_sequence)
			with open(f'{linuxhome_dir}/{project_name}/tmp/{self.idx}_{self.parent_idx}.pdb', 'w') as fout:
				fout.write(pdb+'\n')

		os.system(f"/home/sam/miniconda3/envs/biolib/bin/mkdssp -i {self.idx}_{self.parent_idx}.pdb -o {self.idx}_{self.parent_idx}.ssp")

		#Extract secondary structure and exposure strings
		Seq, Str, Exp = self.ReadSSP(f'{self.idx}_{self.parent_idx}.ssp', exp_threshold, phenotype_level)
		self.ss_structure = Str
		self.exposure = Exp

		#Extract stability
		Stab = self.ReadPDB(f'{self.idx}_{self.parent_idx}.pdb')
		self.stability = Stab

		#Calculate structural complexity.
		self.CalcComplexity()

	def SimilarityToPhenotype(self, target):
		# Currently using Hamming distance but could also implement alignment.
		f = 0
		for p, q in zip(self.ss_structure, target):
			if p == q:	f += 1
		return f

	def AssignFitness(self):

		if fitness_criterion == 'neutral':
			self.fitness = 100

		elif fitness_criterion == 'target_structure':
			self.fitness = self.SimilarityToPhenotype(target_structure)

		elif fitness_criterion == 'complexity_structure':
			self.fitness = self.complexity

		elif fitness_criterion == 'simplicity_structure':
			self.fitness = - self.complexity

	def PrintString(self):
		return f"{self.idx}\t{self.parent_idx}\t{self.nt_sequence}\t{self.aa_sequence}\t{self.ss_structure}\t{self.exposure}\t{self.stability}\t{self.complexity}\t{self.fitness}"
#!/usr/bin/python

# Mat4Pep-matrix_mutator.py
# Jeremy A Horst, 04/04/2010
# updates through 03/09/2018

########################################################################
# INPUT:   a directory of FASTA files with a desired function          #
#          a directory of FASTA files with an undesired function       #
#          a position specific scoring matrix                          #
# GOAL:    tune parameters to maximize score differences between the 2 #
# PROCESS: grid search for gap opening and extension penalty variables #
#          adjust matrix values to maximize                            #
# OUTPUT:  trained matrix, ROC AUC throughout training                 #
#                                                                      #
# requirements:                                                        #
# - python2.5                                                          #
# - FASTA software package                                             #
########################################################################

##################
AA_order = "A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X".split()
matrix_size = len(AA_order)
##################
from sys import argv
from os import mkdir
from os import listdir
from random import shuffle
from random import randrange
from numpy import argmax
from subprocess import call
from subprocess import Popen
from subprocess import PIPE
##################
#make process dir
try: mkdir('tmp')
except: nada=True
##################


######################
# I/O File <> Matrix #
######################
def matrix_file_to_dict(base_substitution_matrix):
	# instantiate dictionary to hold matrix
	matrix_dict = {}
	for AA in AA_order:
		matrix_dict[AA]={}
		for aa in AA_order:
			matrix_dict[AA][aa]=0
	# load in matrix
	for line in open(base_substitution_matrix).readlines():
		if not line.startswith('#') and line.strip():
			if line[0].strip():
				wt_AA = line.split()[0]
				if wt_AA in AA_order:
					for mut_aa in range(matrix_size):
						matrix_dict[wt_AA][AA_order[mut_aa]] = float(line.split()[mut_aa+1])
	return matrix_dict

######################
def write_matrix_from_dict(matrix_dict, matrix_file_name):
	spacer = '   '
	writer = open(matrix_file_name,'w')
	writer.write('    ' + spacer.join(AA_order)+'\n')
	for AA in AA_order:
		printline = AA+' '
		for aa in AA_order:
			value = str(  round(matrix_dict[AA][aa],2)  )[:4]
			if value[0]=='-':  str(  round(matrix_dict[AA][aa],2)  )[:5]
			printline += ' '*(5-len( value )) + value + ' '
		writer.write(printline.strip()+'\n')
	writer.close()
######################


######################
# GAP PENALTY SEARCH #
######################
def train_gap_penalties(matrix_file_name):
	###############################################
	# gop = gap open penalties:       -16, 0, +1  #
	# gep = gap extension penalties:   -8, 0, +1  #
	###############################################
	report_count=0
	first = True
	gap_open_n_extend = [0,0]
	for gap_open in range(-16,0):
		for gap_extension in range(-8,0):
			gap_open      = str(gap_open)
			gap_extension = str(gap_extension)
			# run TSSs
			if   TSS_min_diff:   this_combination = inverter*TSS_minimum_difference(gap_open, gap_extension, matrix_file_name)
			elif TSS_distr_mult: this_combination = inverter*TSS_distribution_mult(gap_open, gap_extension, matrix_file_name)
			elif TSS_distr_add:  this_combination = inverter*TSS_distribution_add(gap_open, gap_extension, matrix_file_name)
			elif TSS_AUC:        this_combination = inverter*TSS_ROC_AUC(gap_open, gap_extension, matrix_file_name)
			if first:
				best_combination = this_combination
				gap_open_n_extend = [gap_open, gap_extension]
				first = False
			if this_combination > best_combination:
				best_combination = this_combination
				gap_open_n_extend = [gap_open, gap_extension]
			report_count+=1
			if gap_open_n_extend[0]:
				print 'gap training... TSS full calculation iterations complete:',report_count,'\t',gap_open_n_extend[0],gap_open_n_extend[1],'\t',gap_open, gap_extension,' '*(4-len(str(gap_extension))),'\tbest_combination:   ',round((best_combination+1)/2,5),'\tthis_combination',round((this_combination+1)/2,5)
			else:
				print 'gap training... TSS full calculation iterations complete:',report_count,'\t',gap_open_n_extend[0],gap_open_n_extend[1],'\t',gap_open, gap_extension
			if best_combination == 1:
				print "Further training is unnecessary; separation is already perfect."
				exit()
	return gap_open_n_extend
######################


#######################################
# TSS STRONG/WEAK COMBINATION OPTIONS #
#######################################
def TSS_minimum_difference(gap_open, gap_extension, matrix_file_name):
	# the goal is to separate the strong from weak binders
	# yet maximizing between worst scoring strong & best scoring weak
	# will skew to this peptides AND
	# AND allow for the inconsistencies of ggsearch
	# so... maximize between the THIRD worst strong & THIRD best weak
	TSS_strong = []
	for high_fasta in high_fastas:
		TSS_ss = TSS([high_fasta], high_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_sw = TSS([high_fasta], high_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_strong += [  TSS_ss - TSS_sw  ]
	TSS_strong.sort(reverse=True)
	third_BIGGEST_strong_TSS = TSS_strong[2]
	TSS_weak = []
	for low_fasta in low_fastas:
		TSS_ws = TSS([low_fasta], low_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_ww = TSS([low_fasta], low_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_weak += [  TSS_ws - TSS_ww  ]
	TSS_weak.sort()
	third_SMALLEST_weak_TSS = TSS_weak[2]
	return (third_SMALLEST_weak_TSS   - third_BIGGEST_strong_TSS) * 1000

#######################################
def TSS_distribution_mult(gap_open, gap_extension, matrix_file_name):
	TSS_ss = TSS(high_fastas, high_lib, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_sw = TSS(high_fastas, high_lib,  low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ws = TSS(low_fastas,  low_lib,  high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ww = TSS(low_fastas,  low_lib,   low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	TSScore = (TSS_ss * TSS_ww) / (TSS_sw * TSS_ws)
	return TSScore

#######################################
def TSS_distribution_add(gap_open, gap_extension, matrix_file_name):
	TSS_ss = TSS(high_fastas, high_lib, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_sw = TSS(high_fastas, high_lib,  low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ws = TSS(low_fastas,  low_lib,  high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ww = TSS(low_fastas,  low_lib,   low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	TSScore = TSS_sw + TSS_ws - TSS_ss - TSS_ww
	return TSScore

#######################################
def TSS_ROC_AUC(gap_open, gap_extension, matrix_file_name):
	TSS_strong = []
	for high_fasta in high_fastas:
		TSS_ss = TSS([high_fasta], high_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_sw = TSS([high_fasta], high_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_strong += [  TSS_ss - TSS_sw  ]
	TSS_weak = []
	for low_fasta in low_fastas:
		TSS_ws = TSS([low_fasta], low_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_ww = TSS([low_fasta], low_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_weak += [  TSS_ws - TSS_ww  ]
	TSS_strong.sort(reverse=True)
	TSS_weak.sort()
	return calc_ROC_AUC(TSS_strong, TSS_weak)
#######################################

def report_TSS_scores(gap_open, gap_extension, matrix_file_name):
	TSS_strong = []
	for high_fasta in high_fastas:
		TSS_ss = TSS([high_fasta], high_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_sw = TSS([high_fasta], high_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_strong += [  TSS_ss - TSS_sw  ]
	TSS_weak = []
	for low_fasta in low_fastas:
		TSS_ws = TSS([low_fasta], low_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_ww = TSS([low_fasta], low_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_weak += [  TSS_ws - TSS_ww  ]
	TSS_strong.sort(reverse=True)
	TSS_weak.sort()
	print '\nstrong:', [x*inverter for x in TSS_strong]
	print '\nweak:',   [x*inverter for x in TSS_weak]
#######################################

def run_ggsearch(gap_open, gap_extension, substitution_matrix, ggsearch_output, fasta, library):
	if method == 'nw':
		command = "%s/bin/ggsearch35_t -C %s -T %s -H -z -1 -d 0 -q -p -f %s -g %s -s %s -O %s %s %s"\
		% (ggsearch_dir, sequence_name_length, threads, gap_open, gap_extension, substitution_matrix, ggsearch_output, fasta, library)
	else:
		command = "%s/bin/ssearch35_t -C %s -T %s -H -z -1 -d 0 -q -p -f %s -g %s -s %s -O %s %s %s"\
		% (ggsearch_dir, sequence_name_length, threads, gap_open, gap_extension, substitution_matrix, ggsearch_output, fasta, library)
	crap = Popen(command.split(),stderr=PIPE,stdout=PIPE).communicate()[0]
#######################################


###################
# TSS CALCULATION #
###################
def TSS(fastasA, setA, fastasB, setB, gap_open, gap_extension, substitution_matrix):
	# TSS_A-B([A]_NA - [B]_NB) = 
	# = (1/ (NA*delta_AB)) sum{1toNA}[ sum{1toNB}[ PSS_ij(1 - delta_ij * delta_AB)]]
	# we add normalization of the second term by the sequence length (not considered before)
	a_names = [entry.split()[0] for entry in open(setA).read()[1:].split('\n>')]
	a_seqs = [''.join(entry.split('\n')[1:]).replace('\n','') for entry in open(setA).read().split('\n>')]
	a_names_seqs = dict(zip(a_names,a_seqs))
	b_names = [entry.split()[0] for entry in open(setB).read()[1:].split('\n>')]
	b_seqs = [''.join(entry.split('\n')[1:]).replace('\n','') for entry in open(setB).read().split('\n>')]
	######################
	# calculate first_term
	# Kronecker delta: check if sets are the same
	delta_AB = 0
	if a_seqs == b_seqs:  delta_AB = 1
	first_term = 1/float(len(a_seqs)*(len(b_seqs)-delta_AB))
	######################
	# for second term iterate through sets
	second_term = 0
	protein_scores = {}
	for protein_a in fastasA:
		protein_a_name = open(setA).read()[1:].split()[0]
		################
		# run 
		ggsearch_output = 'tmp/'+protein_a_name+'.ggsearch_output'
		run_ggsearch(gap_open, gap_extension, substitution_matrix, ggsearch_output, protein_a, setB)
		################
		# parse ggsearch output for each protein in library set
		protein_scores[protein_a_name] = {}
		for line in open(ggsearch_output).read().split('The best scores are:')[1].split('>>>')[0].strip().split('\n')[1:]:
			# collect protein scores
			#sorry, this next line is a parsing hack
			if ')' in line and not ':' in line:
				protein_b_name = line.split()[0]
				score = float(line.split(')')[1].split()[0])
				if not protein_scores[protein_a_name].has_key(protein_b_name):
					protein_scores[protein_a_name][protein_b_name] = score
			elif not line.strip():
				break
	######################
	max_score = max([max(protein_scores[protein_a].values()) for protein_a in protein_scores])
	min_score = min([min(protein_scores[protein_a].values()) for protein_a in protein_scores])
	if method=='sw':  penalty = 0
	else: penalty = calc_penalty(min_score, max_score)
	######################
	for protein_a_name in a_names:
		seq_length = len(a_names_seqs[protein_a_name]) 
		# score protein a vs setB
		for protein_b_name in b_names:
			# Kronecker delta: check if proteins are the same, or even from same protein (uniprot code = 6 alpha-numeric long)
			delta_ab = 0
			if protein_a_name[:unique_sequence_name_length] != protein_b_name[:unique_sequence_name_length]:
				# if no alignment output, assign a bad score
				PairSimScore_ab = penalty
				if protein_scores[protein_a_name].has_key(protein_b_name):
					PairSimScore_ab = protein_scores[protein_a_name][protein_b_name]
				###########################
				# normalize to seq length #
				#   sum to second term    #
	#			second_term += PairSimScore_ab * (1 - delta_ab * delta_AB) / seq_length
				second_term += PairSimScore_ab / seq_length
				###########################
	return first_term * second_term
###################


###################
# TSS CALCULATION #
###################
def calc_ROC_AUC(binders, nonbinders):
	values = [[b,1] for b in binders]
	values += [[n,0] for n in nonbinders]
	total_cases = len(binders)
	total_controls = len(nonbinders)
	# sort (value:answer) tuples by value
	values.sort(reverse=True)
	lenVals = float(len(values))
	# note each change in value & answer as fraction of 0:100
	cases=0
	controls=0
	ROC_exhausting=[]
	for i in values:
		if i[1]: cases += 1
		else:	 controls += 1
		ROC_exhausting+=[[ str(round(100*float(controls)/total_controls,3)) , str(round(100*float(cases)/total_cases,3)) , str(i[0]) ]]
	# ignore changes of only one (value or answer), 
	# catch value before & after change of both.
	## capture the corners
	horiz_vert_change=''
	hold=['0.0','0.0','0.0']
	ROC_terse=[hold]
	for point in ROC_exhausting:
		if point[0][:2]!=hold[0][:2]:
		#if point[0]!=hold[0]:
			if horiz_vert_change=='v':
				ROC_terse+=[hold]
				ROC_terse+=[point]
			horiz_vert_change='h'
		if point[1][:2]!=hold[1][:2]:
		#if point[1]!=hold[1]:
			if horiz_vert_change=='h':
				ROC_terse+=[hold]
				ROC_terse+=[point]
			horiz_vert_change='v'
		hold = point # we keep the last data point before the corner to maintain horizontal & vertical stepwise lines
	ROC_terse+=[point]
	# print AUC
	AUC=0.0
	held_y = 0.0
	for i in ROC_terse:
		x= 1 - float(i[0])/100
		y= float(i[1])/100 - held_y
		AUC += x*y
		held_y = float(i[1])/100
	return 2*AUC-1
def calc_penalty(min_score, max_score):
	# if no alignment was returned, its got to be worse than the worst reported.
	return min_score - ((max_score-min_score)/2)
def directory_to_library(fasta_set,library):
	writer = open(library,'w')
	for f in fasta_set:
		for line in open(f).readlines():
			writer.write(line.strip()+'\n')
	writer.close()
###################


################
# TRAIN MATRIX #
################
def matrix_evolution(base_substitution_matrix, gap_open, gap_extension):
	###################################
	#              CONCEPT            #
	###################################
	# columns & rows mutated en mass to catch large tendancies of amino acid types
	# columns mutated before rows to catch features of the compared target
	#         i.e. the S & W of ss & sw
	# order of individual columns/rows/cells to follow:
	# first follow greedy trajectory by checking each at each step
	# then random local perturbations
	# then monte carlo perturbations
	###################################
	# load matrix file into dictionary
	matrix_dict = matrix_file_to_dict(base_substitution_matrix)
	# write matrix dictionary into file
	matrix_file_name = base_substitution_matrix.split('/')[-1].split('.mat')[0]+'_revised.mat'
	write_matrix_from_dict(matrix_dict, matrix_file_name)
	write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
	###########################
	# calculate base accuracy #
	###########################
	if   TSS_min_diff:   TSS_base = inverter*TSS_minimum_difference(gap_open, gap_extension, matrix_file_name)
	elif TSS_distr_mult: TSS_base = inverter*TSS_distribution_mult(gap_open, gap_extension, matrix_file_name)
	elif TSS_distr_add:  TSS_base = inverter*TSS_distribution_add(gap_open, gap_extension, matrix_file_name)
	elif TSS_AUC:        TSS_base = inverter*TSS_ROC_AUC(gap_open, gap_extension, matrix_file_name)
	TSS_laststep = TSS_base
	
	#########################
	#@# mutate dictionary #@#
	#########################
	
	#########################################
	# start with a deterministic projection #  i.e. GREEDY
	#########################################
	print 'greedy path...'
	if columns_then_cells_OR_cells_only == 0:
		# mutate columns in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 0)
		# mutate rows in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 1)
		# mutate all individual cells in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
		# (this includes local minimization, but if running columns & rows, can do again)
	elif columns_then_cells_OR_cells_only == 1:
		# mutate all individual cells in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
		# (this includes local minimization, but if running columns & rows, can do again)
	
	##########################################
	# proceed with random local minimization #
	##########################################
	if columns_then_cells_OR_cells_only == 0:
		print 'local minimization...'
		# mutate columns in random order
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 0)
		if TSS_AUC and TSS_laststep == 1:
			write_matrix_from_dict(matrix_dict, matrix_file_name)
			print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
			exit()
		# mutate rows in random order
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 1)
		if TSS_AUC and TSS_laststep == 1:
			write_matrix_from_dict(matrix_dict, matrix_file_name)
			print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
			exit()
	
	########################################
	# finish with Monte Carlo minimization #
	########################################
	print 'monte carlo minimization...'
	# repeat with Monte Carlo search
	if columns_then_cells_OR_cells_only == 0:
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 0)
		if TSS_AUC and TSS_laststep == 1:
			write_matrix_from_dict(matrix_dict, matrix_file_name)
			print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
			exit()
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 1)
		if TSS_AUC and TSS_laststep == 1:
			write_matrix_from_dict(matrix_dict, matrix_file_name)
			print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
			exit()
		# mutate all individual cells in random order #
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
		if TSS_AUC and TSS_laststep == 1:
			write_matrix_from_dict(matrix_dict, matrix_file_name)
			print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
			exit()
	elif columns_then_cells_OR_cells_only == 1:
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
		if TSS_AUC and TSS_laststep == 1:
			write_matrix_from_dict(matrix_dict, matrix_file_name)
			print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
			exit()
	return matrix_file_name, TSS_laststep
################
def deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, row_or_column_or_cell):
	# try mutating all columns/rows/cells
	# return to original
	# follow greedy path
	type_crc = ['column','row','cell']
	count_iterations = 0
	if row_or_column_or_cell < 2:
		aa=False
		AAs_done = []
		for i in range(matrix_size):
			# find best AA to improve
			# follow trajectory until end of path
			# find next best AA to improve, repeat
			# go through all AAs, only ONCE (AAs_done)
			increment = 0.04  #this is similar to increments of 1 with min=-19 and max=19
			greedy_TSS        = []
			greedy_trajectory = []
			# find best AA to improve
			for AA in AA_order:
				if not AA in AAs_done:
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					TSS_up,matrix_dict = check_row_column_cell_change(+1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					TSS_down,matrix_dict = check_row_column_cell_change(-1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					try:
						trajectory = 0;  TSS_nextstep = 0
						if TSS_up > TSS_down:  trajectory = +1;  TSSdif = TSS_up
						else:                  trajectory = -1;  TSSdif = TSS_down
						greedy_TSS        += [TSSdif]
						greedy_trajectory += [trajectory]
					except:
						print "PROBLEM: TSS_up, TSS_down, TSSdif, trajectory:",TSS_up, TSS_down, TSSdif, trajectory
				else:
					greedy_TSS        += [0]
					greedy_trajectory += [0]
			if len(AAs_done) < matrix_size:
				# follow trajectory until end of path
				AA           = AA_order[argmax(greedy_TSS)]
				trajectory   = greedy_trajectory[argmax(greedy_TSS)]
				TSS_nextstep = greedy_TSS[argmax(greedy_TSS)]
				if TSS_AUC:
					print AA,'is the best next',type_crc[row_or_column_or_cell],' at step',len(AAs_done)+1,':',round((TSS_nextstep+1)/2,5),'vs',round((TSS_laststep+1)/2,5)
				else:
					print AA,'is the best next',type_crc[row_or_column_or_cell],' at step',len(AAs_done)+1,':',round(TSS_nextstep,5),'vs',round(TSS_laststep,5)
				# still we only follow the trajectory if there is improvement
				if TSS_nextstep > TSS_laststep:
					crap,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, False)
					while TSS_nextstep > TSS_laststep:
						# save
						write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
						TSS_laststep = TSS_nextstep
						# move along trajectory
						TSS_nextstep,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						count_iterations += 1
						if TSS_AUC:
							print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round((TSS_laststep+1)/2,5)
							if TSS_laststep == 1:
								write_matrix_from_dict(matrix_dict, matrix_file_name)
								print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
								exit()
						else:
							print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round(TSS_laststep,5)
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					AAs_done += [AA]
				else:
					# skip the rest; there will be no improvement
					AAs_done += [AA]
					for aaaa in AA_order:
						if not aaaa in AAs_done:
							AAs_done += [aaaa]
	# instance of one amino acid pair at a time
	elif row_or_column_or_cell == 2:
		print 'cell time'
		AA_combos = []
		for AA in AA_order:
			for aa in AA_order:
				if not [aa,AA] in AA_combos:
					AA_combos += [[AA,aa]]
		AAs_done = []
		how_many_we_will_do = len(AA_combos)
		for i in range(how_many_we_will_do):
			if len(AAs_done) < how_many_we_will_do:
				# find best AA combination to improve
				# follow trajectory until end of path
				# find next best AA to improve, repeat
				# go through all AAs, only ONCE (AAs_done)
				increment = 0.04  #this is most similar to -19 to 19
				greedy_TSS        = []
				greedy_trajectory = []
				# find best AA to improve
				for AA,aa in AA_combos:
					if not [AA,aa] in AAs_done and not [aa,AA] in AAs_done:
						matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
						TSS_up,matrix_dict = check_row_column_cell_change(+1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
						TSS_down,matrix_dict = check_row_column_cell_change(-1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
						trajectory = 0;  TSS_nextstep = 0
						if TSS_up > TSS_down:  trajectory = +1;  TSSdif = TSS_up
						else:  trajectory = -1;  TSSdif = TSS_down
						if TSSdif > TSS_laststep:
							if TSS_AUC:
								print AA,aa,(TSS_laststep+1)/2,(TSSdif+1)/2
							else:
								print AA,aa,TSS_laststep,TSSdif
						greedy_TSS        += [TSSdif]
						greedy_trajectory += [trajectory]
					else:
						greedy_TSS        += [0]
						greedy_trajectory += [0]
				if not len(AA_combos) == len(greedy_TSS) or not len(AA_combos) == len(greedy_trajectory):
					print "len(AA_combos) == len(greedy_TSS), len(AA_combos) == len(greedy_trajectory)"
					print len(AA_combos) == len(greedy_TSS), len(AA_combos) == len(greedy_trajectory)
					print "len(AA_combos), len(greedy_TSS), len(AA_combos), len(greedy_trajectory)"
					print len(AA_combos), len(greedy_TSS), len(AA_combos), len(greedy_trajectory)
				# follow trajectory until end of path
				AA,aa        = AA_combos[        argmax(greedy_TSS)]
				trajectory   = greedy_trajectory[argmax(greedy_TSS)]
				TSS_nextstep = greedy_TSS[       argmax(greedy_TSS)]
				if TSS_AUC:
					print AA,aa,'is the best at step:',len(AAs_done)+1,(TSS_laststep+1)/2,(TSS_nextstep+1)/2
				else:
					print AA,aa,'is the best at step:',len(AAs_done)+1,TSS_laststep,TSS_nextstep
				# still we only follow the trajectory if there is improvement
				if TSS_nextstep > TSS_laststep:
					crap,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, False)
					while TSS_nextstep > TSS_laststep:
						# save
						write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
						TSS_laststep = TSS_nextstep
						TSS_nextstep,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						count_iterations += 1
						if TSS_AUC:
							print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round((TSS_laststep+1)/2,5)
							if TSS_laststep == 1:
								write_matrix_from_dict(matrix_dict, matrix_file_name)
								print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
								exit()
						else:
							print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round(TSS_laststep,5)
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					AAs_done += [[AA,aa]]
					if not AA==aa:  AAs_done += [[aa,AA]]
				else:
					# skip the rest; there will be no improvement
					AAs_done += [[AA,aa]]
					if not AA==aa:  AAs_done += [[aa,AA]]
					for AAAaaa in AA_combos:
						if not AAAaaa in AAs_done:
							AAs_done += [AAAaaa]
	return matrix_dict, TSS_laststep
################
def check_row_column_cell_change(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, run_TSS):
	# columns
	if row_or_column_or_cell == 0:
		for aa in AA_order:
			matrix_dict[aa][AA] = check_to_maxmin(matrix_dict[aa][AA], increment)
	# row
	elif row_or_column_or_cell == 1:
		for aa in AA_order:
			matrix_dict[AA][aa] = check_to_maxmin(matrix_dict[AA][aa], increment)
	# cells
	elif row_or_column_or_cell == 2:
		matrix_dict[AA][aa] = check_to_maxmin(matrix_dict[AA][aa], increment)
		if not AA==aa:
			matrix_dict[aa][AA] = check_to_maxmin(matrix_dict[aa][AA], increment)
	# calc
	write_matrix_from_dict(matrix_dict, matrix_file_name)
	if run_TSS:
		if   TSS_min_diff:   TSS = inverter*TSS_minimum_difference(gap_open, gap_extension, matrix_file_name)
		elif TSS_distr_mult: TSS = inverter*TSS_distribution_mult(gap_open, gap_extension, matrix_file_name)
		elif TSS_distr_add:  TSS = inverter*TSS_distribution_add(gap_open, gap_extension, matrix_file_name)
		elif TSS_AUC:        TSS = inverter*TSS_ROC_AUC(gap_open, gap_extension, matrix_file_name)
	else: TSS = 0
	return TSS, matrix_dict
################
def check_to_maxmin(valueA, increment):
	out = 0
	if valueA+increment > value_max:
		out = value_max
	elif valueA+increment < value_min:
		out = value_min
	else:
		out = valueA+increment
	return out
################
def mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, row_or_column_or_cell):
	increment = 0.04  #this is similar to increments of 1 with min=-19 and max=19
	count_iterations = 0
	# mutate rows in random order
	if row_or_column_or_cell < 2:
		row_order  = []
		row_order += AA_order
		shuffle(row_order)
		for AA in row_order:
			matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, False)
	elif row_or_column_or_cell == 2:
		AA_matrix = []
		for AA in AA_order:
			for aa in AA_order:
				if not [aa,AA] in AA_matrix:
					AA_matrix += [[AA,aa]]
		shuffle(AA_matrix)
		for (AA,aa) in AA_matrix:
			matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, aa)
	return matrix_dict, TSS_laststep
################
def mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, row_or_column_or_cell):
	count_iterations = 0
	###############################
	# mutate rows in random order #
	###############################
	if row_or_column_or_cell < 2:
		row_order  = []
		row_order += AA_order
		shuffle(row_order)
		for AA in row_order:
		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
		#$$$$$  Bounce out to random distance  $$$$$#
		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
			for i in range(5):
				# increment should not be zero, but pos & neg are vital
				increment = randrange(43)/10.0	#this is most similar to -19 to 19 / 10
				posneg = randrange(2)
				if posneg: increment = -1*increment
				# run...
				matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, False)
	elif row_or_column_or_cell == 2:
		AA_matrix = []
		for AA in AA_order:
			for aa in AA_order:
				if not [aa,AA] in AA_matrix:
					AA_matrix += [[AA,aa]]
		shuffle(AA_matrix)
		for (AA,aa) in AA_matrix:
			for i in range(5):
				increment = randrange(43)/10.0	#this is most similar to -19 to 19 / 10
				posneg = randrange(2)
				if posneg: increment = -1*increment
				# run...
				matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, aa)
	return matrix_dict, TSS_laststep
################
def increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, aa):
	#@# first, get trajectory
	type_crc = ['column','row','cell']
	###################
	# save dictionary #
	###################
	matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
	##################
	# try adding one #
	##################
	TSS_up,matrix_dict = check_row_column_cell_change(+1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
	matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
	#######################
	# try subtracting one #
	#######################
	TSS_down,matrix_dict = check_row_column_cell_change(-1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
	matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
	###########
	# compare #
	###########
	trajectory = 0;  TSS_nextstep = 0
	if   TSS_up   > TSS_down:  trajectory = +1;  TSS_nextstep = TSS_up
	else:  trajectory = -1;  TSS_nextstep = TSS_down
	count_iterations += 1
	if TSS_AUC:
		print AA,aa,count_iterations,'iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round((TSS_laststep+1)/2,5),'Next_step:',round((TSS_nextstep+1)/2,5)
		if TSS_laststep == 1:
			write_matrix_from_dict(matrix_dict, matrix_file_name)
			print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
			exit()
	else:
		print AA,aa,count_iterations,'iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round(TSS_laststep,5),'Next_step:',round(TSS_nextstep,5)
	###########################################################
	#@# now, move along trajectory unil no more improvement #@#
	###########################################################
	if TSS_nextstep > TSS_laststep:
		# push back to improved trajectory; the next line will proceed... n+1+1-1=n+1
		matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
		crap,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, False)
		while TSS_nextstep > TSS_laststep:
			# save
			write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
			TSS_laststep = TSS_nextstep
			TSS_nextstep,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
			count_iterations += 1
			if TSS_AUC:
				print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round((TSS_laststep+1)/2,5)
				if TSS_laststep == 1:
					write_matrix_from_dict(matrix_dict, matrix_file_name)
					print "Further training is unnecessary; separation is already perfect. Trained matrix:",matrix_file_name
					exit()
			else:
				print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',round(TSS_laststep,5)
		# revert to most recent BTTR form
		matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
	return matrix_dict, TSS_laststep, count_iterations
################

		#########
		# START #
		#########

if __name__=='__main__':
	try:
		######################
		# prepare input sets #
		######################
		# input directories
		high_dir = argv[1]+'/'
		low_dir = argv[2]+'/'
		# input fasta lists
		high_fastas=[]
		for f in listdir(high_dir):
			if f.endswith('.fasta'): high_fastas+=[high_dir+f]
		low_fastas=[]
		for f in listdir(low_dir):
			if f.endswith('.fasta'): low_fastas+=[low_dir+f]
		# make searchable libraries
		high_lib = "lib_high.fasta"
		low_lib  = "lib_low.fasta"
		directory_to_library(high_fastas,high_lib)
		directory_to_library(low_fastas,low_lib)
		
		###########
		# options #
		###########
		columns_then_cells_OR_cells_only = 0
		if '-cells' in argv:  columns_then_cells_OR_cells_only = 1
		substitution_matrix = 'pam250.mat'
		if '-matrix' in argv:  substitution_matrix = argv[argv.index('-matrix')+1]
		ggsearch_dir = './fasta-35.4.11/'
		if '-fasta_dir' in argv:  ggsearch_dir = argv[argv.index('-fasta_dir')+1]
		method = 'nw'
		if '-sw' in argv: method = 'sw'
		threads = '2'
		if '-threads' in argv:  threads = argv[argv.index('-threads')+1]
		sequence_name_length = 8
		if '-namelen' in argv:  sequence_name_length = int(argv[argv.index('-namelen')+1])
		unique_sequence_name_length = 6
		if '-uniqname' in argv:  unique_sequence_name_length = int(argv[argv.index('-uniqname')+1])
		allow_inversion = False
		if '-invert' in argv: allow_inversion=True
		
		TSS_distr_add  = True; TSS_min_diff = False; TSS_distr_mult = False; TSS_AUC = False
		if   '-multiply' in argv: TSS_distr_add  = False; TSS_min_diff = False; TSS_distr_mult = True; TSS_AUC = False
		elif '-maximize_diff' in argv: TSS_distr_add  = False; TSS_min_diff = True; TSS_distr_mult = False; TSS_AUC = False
		elif '-ROC_AUC' in argv: TSS_distr_add  = False; TSS_min_diff = False; TSS_distr_mult = False; TSS_AUC = True
		
		gaps_only = False
		if '-gaps_only' in argv: gaps_only = True  #train gap open values ONLY
		gap_open=False; gap_extension=False
		if '-gop_gep' in argv: # set gap open & extend values
			gap_open = int(argv[argv.index('-gop_gep')+1])
			gap_extension = int(argv[argv.index('-gop_gep')+2])
			if gap_open>0 or gap_extension>0:
				print "gap values should be integers < 0."
				exit()
			
		###############################################################
		# the above names & lists will be used throughout the process #
		###############################################################
	except:
		print "\nUsage: matrix_mutator.py <desired_peptides/> <unwanted_peptides/>"
		print "Options:"
		print "\t-matrix <matrix>\tBase substitution matrix (default=pam250)."
		print "\t-fasta_dir <dir>\tPath to the FASTA suite.\n"
		
		print "\t-cells\t\t\tMutate cell by cell, trajectory, & then monte carlo."
		print "\t\t\t\tDefault is trajectory by column,row,cell; then monte carlo by column,row.\n"
		print "\t-multiply\t\tDrive matrix training with ss*ww/(sw*ws) [s=strong=functional, w=weak=no_fxn]."
		print "\t-maximize_diff\t\tMaximize separation of 3rd highest scoring weak & 3rd lowest scoring strong peptide."
		print "\t-ROC_AUC\t\tTrain to Area Under the Reciever Operator Characteristic Curve (default=TSS)."
		print "\t\t\t\tDefault matrix mutation driving force is ss+ww-sw-ws.\n"
		
		print "\t-sw\t\t\tApply Smith-Waterman scoring (default=Needleman-Wunsch)."
		print "\t-invert\t\t\tAllow inversion: multiply output by -1. Useful if matrices separate peptides in wrong direction, but validity is arguable."
		print "\t-gaps_only\t\tOnly do grid search for gap values, do not alter matrix."
		print "\t-gop_gep # #\t\tInput gap values (integers < 0).\n"
		
		print "\t-namelen #\t\tCatches sequence names, needed due to ggsearch (default=8)."
		print "\t-uniqname #\t\tAvoids comparing peptides from the same protein, assumes UniProt coding (default=6)."
		print "\t-threads #\t\tUse more than 1 processor (default=2)."
	
	# get & set min & max values. Will be used to contain run-away training
	string_values = [line.strip().split()[1:] for line in open(substitution_matrix).readlines() if not line[0]=='#' and line[0].strip()]
	min_value = min([float(x) for line in string_values for x in line])
	max_value = max([float(x) for line in string_values for x in line])
	value_min = max_value - (max_value - min_value)*1.5
	value_max = (max_value - min_value)*1.5 + min_value
	####################################################################
	
	####################################################################
	# Check whether to invert by multiplying output by -1.
	inverter = 1
	if allow_inversion:
		score = TSS_ROC_AUC(-8, -4, substitution_matrix)
		if score < 0:  inverter = -1
		print ("", "Do not invert", "Invert")[inverter], substitution_matrix, "to achieve best separation.",(score*inverter+1)/2
		if score*inverter == 1:
			print "Training this matrix is unnecessary; separation is already perfect."
			report_TSS_scores(-8, -4, substitution_matrix)
			exit()
	else:
		# set gap penalties
		if not gap_open or not gap_extension:
			gap_open, gap_extension = train_gap_penalties(substitution_matrix)
		print substitution_matrix,'  gap open penalty:',gap_open,  'gap extend penalty:',gap_extension
		
		#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
		# mutate each matrix to maximize (TSS_strong-strong - TSS_strong-weak) #
		#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
		if not gaps_only:
			substitution_matrix, TSS_laststep = matrix_evolution(substitution_matrix,   gap_open,    gap_extension)
			if TSS_AUC:  print 'Trained matrix acheives TSS difference of:',(TSS_laststep+1)/2
			else:  print 'Trained matrix acheives TSS difference of:',TSS_laststep
			print 'Trained matrix file:', substitution_matrix
			print substitution_matrix,'  gap open penalty:',gap_open,  'gap extend penalty:',gap_extension, 'multiplier:',inverter

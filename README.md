# Mat4Pep
# Pellitrix
Machine learning tool for functional peptide sets

by Jeremy Horst
2010-05

###############################
run the pellitrix program with:

./matrix_scorer.peptide.py <sequence.fasta> <scoring matrix> [gop] [gep]

use the script "matrix_scorer.peptide.py" to score an entire sequence.

use the script "matrix_scorer.protein.py" to score each residue in the sequence.
################################################################################


For detailed explanation of the underlying science, see:

pelltrix/refs/Pellitrix_manuscript.pdf

For source of the binders, see:

pelltrix/refs/SiqueiraOppenheim2009.pdf


##########################################################################################
This package requires python2.5 (or a later version) and fasta33 to be properly installed.

To get fasta33 working, untar the fasta_pkg.tgz with the following command: 

tar -zxf fasta_prg.tgz

Read the README to install. Basically:

cd fasta-35.4.11/src/

make -f ../make/Makefile.linux32
################################


for questions email: jahorst@gmail.com

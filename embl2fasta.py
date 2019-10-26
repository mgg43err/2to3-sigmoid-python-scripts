import sys
from Bio import SeqIO
SeqIO.convert(sys.stdin, "embl", sys.stdout, "fasta")

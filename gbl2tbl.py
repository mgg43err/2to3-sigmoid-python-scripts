#pip install biopython
#pip install geneparse
import Bio
print(Bio.__version__)

import sys
from Bio import SeqIO
SeqIO.convert(sys.stdin, "embl", sys.stdout, "fasta")

import os, sys, platform
import os.path
import argparse
from Bio import SeqIO
from geneparse import *
def createParser():
  parser = argparse.ArgumentParser(
      prog='Genbank to .tbl converter',
      usage='\n%(prog)s <input_file> [options]',
      description = '''This script allows to convert genbank file into .tbl file format
      Requires Biopython 1.64 (or newer). Use '>' (stdout) to create a .tbl file''',
      epilog='''(c) Aliaksandr Damienikan, 2015.'''
  )
  parser.add_argument('input_file',
                      help='''path to input Genbank file.''')
  parser.add_argument('-f', '--fasta',
                      default=False,
                      action='store_const',
                      const=True,
                      help='''creates fasta from genbank file.''')
  parser.add_argument('-p', '--prefix',
                      default=False,
                      help='''sequencing centre prefix.''')
  parser.add_argument('-t', '--translation',
                      default=False,
                      action='store_const',
                      const=True,
                      help='''adds translation qualifier to CDS features in .tbl''')
  parser.add_argument('-v','--version', action='version', version='%(prog)s 1.3 (October 28, 2017)')
  return parser
args = createParser()
enter = args.parse_args()
arguments = sys.argv[1:0]
cwd = os.path.abspath(os.path.dirname(sys.argv[0])) # or __file__ in hooks ?
genbank_file = SeqIO.parse(enter.input_file, 'genbank')
if enter.fasta != False:
    if platform.system() == 'Windows':
        fasta_name = enter.input_file[enter.input_file.rfind('\\'):enter.input_file.rfind('.')]
        fasta_file = open('%s\\%s.fasta' % (cwd, fasta_name), 'w')
    else:
        fasta_name = enter.input_file[enter.input_file.rfind('/'):enter.input_file.rfind('.')]
        fasta_file = open('%s/%s.fasta' % (cwd, fasta_name), 'w')
    SeqIO.write(genbank_file, fasta_file, 'fasta')
    genbank_file.close()
    fasta_file.close()
genbank_file = GenbankFile(path=enter.input_file)
genbank_file.parse().count = -1
genbank_file = SeqIO.parse(enter.input_file, 'genbank')
for record in genbank_file:
  print('>Feature %s' % record.id)
  count += 1
  genbank = GenbankFile(path=enter.input_file)
  genbank.parse()
  feature_type = str(feature.type)
  qualifiers = ''
  for pair in feature.qualifiers:
    for key, value in pair.items():
      if enter.translation == True and key == 'translation':
        qualifiers += '\n\t\t\t%s\t%s' % (key, value)
        if key == 'protein_id':
          if enter.prefix != False:
            for find in feature.qualifiers:
              for adkey, advalue in find.items():
                if adkey == 'locus_tag':
                  loci = advalue
                  break
                  qualifiers += '\n\t\t\t%s\t%s' % (key, value)
                  if key != 'protein_id' and key !='translation':
                    qualifiers += '\n\t\t\t%s\t%s' % (key, value)
                    print('%s\t%s\t%s\t%s' % (start, end, feature_type, qualifiers))
                    genbank_file.close()

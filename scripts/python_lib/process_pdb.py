import os
import re
import numpy
import sys
import click

# This function will parse a given PDB file and retain only the ATOM records and eliminate all the HETATM records
# Also, in case of alternate locations like ASER and BSER, it will retain the 'A' records only.

@click.command()
@click.option('--input_dir', required=True, type=str, help='Name of the \
	input directory containing the raw pdb files')
@click.option('--output_dir', required=True, type=str, help='Name of the directory \
	to which the processed pdb files will be written to')

def main (input_dir, output_dir):
	last_char = input_dir[-1]
	if last_char != '/':
		input_dir = input_dir + '/'
	last_char = output_dir + '/'
	if last_char != '/':
		output_dir = output_dir + '/'

	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)
	pdbfiles = os.listdir(input_dir)
	for pdbfile in pdbfiles:
		if not re.match('.*?\.pdb$', pdbfile):
			continue
		else:
			pdbfile_fq = input_dir  + pdbfile
			outfile_fq = output_dir + pdbfile
			fh1 = open(pdbfile_fq, 'r')
			fh2 = open(outfile_fq, 'w')
			print("Processing file ", pdbfile)
			for line in fh1:
				if re.match('^ATOM', line):
					alt_loc = line[16] # alternate location if any
					alt_loc = alt_loc.replace(' ','')
					if alt_loc == 'A': # we will consider by default the 'A' location
						fh2.write(line)
					elif alt_loc == '':
						fh2.write(line)
					else:
						continue						
				elif re.match('^TER', line):
					fh2.write(line)
				else:
					continue
			fh1.close()
			fh2.close()
			print("Done!\n")


if __name__ == '__main__':
	main()

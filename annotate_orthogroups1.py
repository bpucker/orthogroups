### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """ python annotate_orthogroups1.py
							--out <OUTPUT_FILE>
							--og <ORTHOGROUP_FILE>
							--in <INPUT_FASTA_FILE_FOLDER>
							feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""

import os, sys, glob

# --- end of imports --- #


def load_multiple_fasta_file( fasta_file ):
	"""!@brief load content of multiple fasta file """
	
	content = {}
	with open( fasta_file, "r" ) as f:
		header = f.readline().strip()[1:]
		if " " in header:
			header = header.split(' ')[0]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				content.update( { header: "".join( seq ).upper() } )
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		content.update( { header: "".join( seq ).upper() } )
	return content


def load_orthogroups( orthofinder_result_file ):
	"""! @brief load orthogroups into nested dictionary """
	
	orthogroups = {}
	with open( orthofinder_result_file, "r" ) as f:
		samples = f.readline().strip().replace( ".pep", "" ).split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.split('\t')
			parts[-1] = parts[-1].strip()
			group = []
			for idx, part in enumerate( parts[1:] ):
				if len( part ) > 1:
					if "," in part:
						subparts = part.split(', ')
					else:
						subparts = [ part ]
				else:
					subparts = []
				for each in subparts:
					group.append( each )
			if len( group ) > 0:
				orthogroups.update( { parts[0]: group } )
			line = f.readline()
	return orthogroups


def main( arguments ):
	
	fasta_file_folder = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	orthogroup_file = arguments[ arguments.index('--og')+1 ]
	
	# --- get IDs per orthogroup --- #
	ogs = load_orthogroups( orthogroup_file )	#dictionary: og name as key and genes in list as value
	
	# --- load all sequences --- #
	fasta_filenames = glob.glob( fasta_file_folder + "*.fasta" ) + glob.glob( fasta_file_folder + "*.pep.fa" )
	sequences = {}
	for filename in fasta_filenames:
		sequences.update( load_multiple_fasta_file( filename ) )
	
	# --- write sequences into output file --- #
	with open( output_file, "w" ) as out:
		for og in ogs.keys():
			members = ogs[ og ]	#list of sequence IDs
			for idx, member in enumerate( members ):
				try:
					out.write( '>' + og + "_%_" + member + "\n" + sequences[ member ] + "\n" )
				except KeyError:
					print member


if '--out' in sys.argv and '--in' in sys.argv and '--og' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

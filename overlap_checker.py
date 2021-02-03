### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """	python overlap_checker.py
							--ref <REFERENCE_ORTHOGROPU>
							--ogs <ORTHOGROUP_FILES>
							--out <OUTPUT_FOLDER>
					"""

import re, math, os, sys, glob
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import DataFrame

# --- end of imports --- #


def load_orthogroups( orthofinder_result_file ):
	"""! @brief load orthogroups into nested dictionary """
	
	orthogroups = {}
	orthogroup_sizes = {}
	with open( orthofinder_result_file, "r" ) as f:
		samples = f.readline().strip().replace( ".pep", "" ).split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.split('\t')
			parts[-1] = parts[-1].strip()
			group = {}
			group_size = 0
			for idx, part in enumerate( parts[1:] ):
				if len( part ) > 1:
					subparts = part.split(', ')
				else:
					subparts = []
				group.update( { samples[ idx ]: subparts } )
				group_size += len( subparts )
			if len( group.keys() ) > 0:
				orthogroups.update( { parts[0]: group } )
				orthogroup_sizes.update( { parts[0]: group_size } )
			line = f.readline()
	return orthogroups, orthogroup_sizes


def load_og_per_gene( og_size_file ):
	"""! @brief load OG per gene """
	
	og_per_gene = {}
	with open( og_size_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			og_per_gene.update( { parts[0]: int( parts[1] ) } )
			line = f.readline()
	return og_per_gene


def load_genes_per_og( orthofinder_result_file1 ):
	"""! @brief load all genes per og """
	
	ogs = []
	names = []
	repr_specs = []
	with open( orthofinder_result_file1, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			genes = []
			spec_counter = 0
			for part in parts[1:]:
				if len( part.strip() ) > 1:
					spec_counter += 1
					if "," in part:
						genes.append( part.split(',') )
					else:
						genes.append( [ part ] )
			genes = [ x for block in genes for x in block ]
			fin_genes = {}
			for gene in genes:
				fin_genes.update( { gene: None } )
			ogs.append( fin_genes )
			names.append( parts[0] )
			repr_specs.append( spec_counter )
			line = f.readline()
	return ogs, names, repr_specs


def calculate_recovery( string ):
	"""! @brief calculate recovered percentage """
	
	if "," in string:
		values = map( int, string.split(',') )
	else:
		values[ int( values ) ]
	if sum( values ):
		return max( values ) / float( sum( values ) )
	else:
		return 0


def construct_figure( data1, data2, data3, figfile ):
	"""! @brief construct a summary figure to show the orthogropu recovery between analyses """
	
	#gene, pigment_state, spec, value
	# datamatrix = [ 	[ "CHS", "A", "x1", 10 ],
								# [ "CHS", "A", "x2", 8 ],
								# [ "CHS", "B", "x3", 5 ],
								# [ "CHS", "B", "x4", 3 ],
								# [ "CHI", "A", "x1", 5 ],
								# [ "CHI", "A", "x2", 7 ]
							# ]
	
	datamatrix = []
	for og in data1:
		datamatrix.append( [ "OG2", og ] )
	for og in data2:
		datamatrix.append( [ "OG3", og ] )
	for og in data3:
		datamatrix.append( [ "OG4", og ] )
	
	df = DataFrame( datamatrix, columns=[ "name", "overlap" ])
					
	# -- generate figure --- #
	fig, ax = plt.subplots( )	#, sharex=True
	
	x = sns.violinplot( 	x="lineage",
									y="gene expression",
									data = df,	#DataFrame: column=variable, row=observation
									palette = [ "dodgerblue", "magenta", "black" ],
									scale = "count",
									cut = 0
								)
	
	ax.text( 0, 1, "n=" + str( len( data1 ) ), va="bottom", ha="center" )
	ax.text( 1, 1, "n=" + str( len( data2 ) ), va="bottom", ha="center" )
	ax.text( 2, 1, "n=" + str( len( data3 ) ), va="bottom", ha="center" )
	
	ax.set_ylabel( "percentage of OG1 recovery", fontsize=universal_fontsize )
	ax.set_xlabel( "dataset", fontsize=universal_fontsize )
	
	ax.tick_params(axis='x', which='major', labelsize=universal_fontsize, rotation=90 )
	ax.tick_params(axis='x', which='minor', labelsize=universal_fontsize, rotation=90 )
	
	ax.tick_params(axis='y', which='major', labelsize=universal_fontsize )
	ax.tick_params(axis='y', which='minor', labelsize=universal_fontsize )
	
	plt.subplots_adjust( left=0.15, right=0.999, top=0.95, bottom=0.15, hspace=0.03 )
	
	fig.savefig( figure_file, dpi=300 )
	try:
		fig.savefig( figfile.replace(".pdf", ".png"), dpi=300 )
	except:
		pass
	plt.close("all")


def main( arguments ):
	"""! @brief run everything """
	
	ref_orthofinder_result_file = arguments[ arguments.index('--ref')+1 ]
	orthofinder_result_files = arguments[ arguments.index('--ogs')+1 ]
	if "," in orthofinder_result_files:
		orthofinder_result_files = orthofinder_result_files.split(',')
	else:
		orthofinder_result_files = [ orthofinder_result_files ]
	sample_names = []
	for each in orthofinder_result_files:
		sample_names.append( each.split('/')[-1].split('.')[0] )
	
	output_folder = arguments[ arguments.index('--out')+1 ]

	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )

	# ---- load data --- #
	ref_ogs, ref_names, ref_repr_specs = load_genes_per_og( ref_orthofinder_result_file )
	
	collected_data = {}
	for k, sample in enumerate( sample_names ):
		ogs, names, repr_specs= load_genes_per_og( orthofinder_result_files[ k ] )
		collected_data.update( { sample: { 'ogs': ogs, 'names': names, 'repr_specs': repr_specs } } )


	# --- check overlap between orthogroups --- #
	overlap_table = output_folder + "overlap_table.txt"
	if not os.path.isfile( overlap_table ):
		with open( overlap_table, "w", 0 ) as out:
			out.write( "RefOG\tOG_size\t" + "\t".join( sample_names ) + "\n" )
			for idx, current_og in enumerate( ref_ogs ):
				new_line = [ ref_names[ idx ], str( len( current_og.keys() ) ) ]
				for s, sample in enumerate( sample_names ):
					matches = []
					for og in collected_data[ sample ]["ogs"]:
						match = len( list( set( og ).intersection( current_og ) ) )
						if match > 0:
							matches.append( match )
					new_line.append( str( ",".join( map( str, matches ) ) ) )
				out.write( "\t".join( new_line )  + "\n" )	#one line per orthogroup

	# --- generate figure --- #
	og2 = []
	og3 = []
	og4 = []

	with open( overlap_table, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.split('\t')
			parts[-1] = parts[-1].strip()	#careful way of stripping in case of empty fields
			
			# --- OG2 analysis --- #
			og2.append( calculate_recovery( parts[1] ) )
			
			# --- OG3 analysis --- #
			og3.append( calculate_recovery( parts[2] ) )
			
			# --- OG4 analysis --- #
			og4.append( calculate_recovery( parts[3] ) )
			
			line = f.readline()

	figfile = output_folder + "orthogroup_overlap.pdf"
	construct_figure( og2, og3, og4, figfile )


if '--ref' in sys.argv and '--ogs' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

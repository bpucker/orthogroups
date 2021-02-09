### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###

__usage__ = """
					python enrichment_check.py
					--in <INPUT_FILE>
					--out <OUTPUT_FOLDER>
					
					optional:
					--anno <ANNOTATION_FILE>
					"""


import re, math, os, sys, glob
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

# --- end of imports --- #


def load_orthogroups( orthofinder_result_file, pigment_states ):
	"""! @brief load orthogroups into nested dictionary """
	
	orthogroups = {}
	orthogroup_sizes = {}
	A_sizes = {}
	B_sizes = {}
	with open( orthofinder_result_file, "r" ) as f:
		samples = f.readline().strip().replace( ".pep", "" ).split('\t')[1:]	#species IDs
		line = f.readline()
		while line:
			parts = line.split('\t')
			parts[-1] = parts[-1].strip()	#remove line break from last field (strip() on line would remove empty columns)
			group = {}
			group_size = 0
			A_size = 0
			B_size = 0
			for idx, part in enumerate( parts[1:] ):
				if len( part ) > 1:
					if "," in part:
						subparts = part.split(', ')
					else:
						subparts = [ part.strip() ]
				else:
					subparts = []
				group.update( { samples[ idx ]: subparts } )
				group_size += len( subparts )
				try:
					pigment = pigment_states[ samples[ idx ] ][0]	#this could also be used to look at individual lineages
					if pigment == "A":
						A_size += len( subparts )
					elif pigment == "B":
						B_size += len( subparts )
				except KeyError:
					pass
			if len( group.keys() ) > 0:
				orthogroups.update( { parts[0]: group } )
				orthogroup_sizes.update( { parts[0]: group_size } )
				A_sizes.update( { parts[0]: A_size } )
				B_sizes.update( { parts[0]: B_size } )
			line = f.readline()
	return orthogroups, orthogroup_sizes, A_sizes, B_sizes


def load_pigment_states( taxon_file ):
	"""! @brief load pigment states of all taxa """
	
	pigment_state = {}
	with open( taxon_file, "r" ) as f:
		f.readline()	#remove header line
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			pigment_state.update( { parts[1]: parts[3] } )
			line = f.readline()
	return pigment_state


def summarize_orthogroup_sizes( A_sizes, B_sizes, figfile ):
	"""! @brief analyze group sizes """
	
	A = A_sizes.values()
	B = B_sizes.values()
	
	fig, ax = plt.subplots()
	
	Ay_values = []
	x_values = []
	By_values = []
	for i in range( 50 ):
		x_values.append( i+1 )
		Ay_values.append( A.count( i+1 ) )
		By_values.append( B.count( i+1 ) )
	
	ax.plot( x_values, Ay_values, linestyle="--", marker="o", color="blue", label="A" )
	ax.plot( x_values, By_values, linestyle="--", marker="o", color="red", label="B" )
	
	ax.set_xlim( 0, 50 )
	
	ax.set_xlabel( "orthogroup size" )
	ax.set_xlabel( "number of orthogroups" )
	
	fig.savefig( figfile, dpi=300 )
	plt.close("all")
	
	#return np.mean( A ), np.mean( B )	
	return np.median( Ay_values ), np.median( By_values )	


def compare_group_sizes_across_pigment_states( A_sizes, B_sizes, fig_file, A_avg, B_avg, l2fc_cutoff ):
	"""! @brief compare contribution to orthogroup size between pigment states """
	
	# --- collect all required data --- #
	data = []
	x_values = []
	y_values = []
	up_x_values = []
	up_y_values = []
	down_x_values = []
	down_y_values = []
	for key in A_sizes.keys():
		A = ( A_sizes[ key ] + 1.0 ) / A_avg
		B = ( B_sizes[ key ] + 1.0 ) / B_avg
		m = math.log( ( A*B ), 10 )
		a = math.log( ( B / A ), 2 )
		data.append( { 'ID': key, 'm': m, 'a': a, 'A': A, 'B': B } )
		
		if a > l2fc_cutoff:
			up_x_values.append( m )
			up_y_values.append( a )
		elif a < -1 * l2fc_cutoff:
			down_x_values.append( m )
			down_y_values.append( a )
		else:
			x_values.append( m )
			y_values.append( a )
		
	
	# --- construct figure --- #
	fig, ax = plt.subplots()
	
	ax.scatter( x_values, y_values, color="grey", s=1, marker="." )
	ax.scatter( up_x_values, up_y_values, color="red", s=1, marker="." )
	ax.scatter( down_x_values, down_y_values, color="blue", s=1, marker="." )
	
	ax.set_xlabel( "relative orthogroup size [log10( A*B )]" )
	ax.set_ylabel( "orthogroup size change in B [log2( B / A )]" )
	
	fig.savefig( fig_file, dpi=300 )
	plt.close( "all" )
	
	return data


def load_annotation( anno_file ):
	"""! @brief load annotation from given file """
	
	anno = {}
	with open( anno_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno.update( { parts[0]: "; ".join( parts[1:] ) } )
			line = f.readline()
	return anno


def main( arguments ):
	"""! @brief run everything """
	
	orthofinder_result_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	taxon_file = arguments[ arguments.index('--taxon')+1 ]
	if '--anno' in arguments:
		anno_file = arguments[ arguments.index('--anno')+1 ]
		anno = load_annotation( anno_file )
	else:
		anno = {}
	
	l2fc_cutoff = 2
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	pigment_states = load_pigment_states( taxon_file )
	orthogroups, og_sizes, A_sizes, B_sizes = load_orthogroups( orthofinder_result_file, pigment_states )
	
	# --- analyze size of orthogroups --- #
	fig_file1 = output_folder + "pigment_state_group_sizes.pdf"
	A_avg, B_avg = summarize_orthogroup_sizes( A_sizes, B_sizes, fig_file1 )
	print "A median: " + str( A_avg )
	print "B median: " + str( B_avg )
	
	# --- compare orthogroup contribution between pigment states --- #
	fig_file2 = output_folder + "pigment_state_group_size_comparison.pdf"
	data = compare_group_sizes_across_pigment_states( A_sizes, B_sizes, fig_file2, A_avg, B_avg, l2fc_cutoff )
	
	# --- data output in table --- #
	output_summary_table = output_folder + "summary.txt"
	with open( output_summary_table, "w" ) as out:
		out.write( "ID\tlog10(OrthoGroupSize)\tlog2(OrthoGroupEnrichment)\tA\tB\tAnno\n" )
		data = sorted( data, key=itemgetter('a') )[::-1]
		for each in data:
			new_line = [ each['ID'], each['m'], each['a'], each['A'], each['B'] ]
			try:
				new_line.append( anno[ each['ID'] ] )
			except KeyError:
				new_line.append( "." )
			out.write( "\t".join( map( str, new_line ) ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv and "--taxon" in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

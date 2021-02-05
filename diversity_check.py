### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """ python diversity_check.py
							--out <OUTPUT_FOLDER>
							--in <INPUT_OG_FILE>
							feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""

import re, os, sys, math
from operator import itemgetter

# --- end of imports --- #

def calculate_shannon_index( values ):
	"""! @brief calculate shannon index based on given value distribution """
	
	values = [ i for i in values if i != 0 ]
	og_size = float( sum( values ) )
	h = []
	for val in values:
		p = val / og_size
		h.append( p * math.log( p, 2 ) )
	shannon = -1 * sum( h )
	eveness = shannon / ( -1 * math.log( 1.0 / len( values ) ) )
	return shannon, eveness


def construct_figure( fig_file, xvalues, yvalues, yname ):
	"""! @brief plot the orthogroup sizes against the shannon index """
	
	import matplotlib.pyplot as plt
	fig, ax = plt.subplots()
	
	ax.scatter( xvalues, yvalues, s=1, color="black" )
	
	ax.set_xlabel( "log10( orthogroup size )" )
	ax.set_ylabel( yname )
	
	fig.savefig( fig_file, dpi=300 )


def main( arguments ):
	""""! @brief run everything """
	
	output_folder = arguments[ arguments.index('--out')+1 ]
	input_file = arguments[ arguments.index('--in')+1 ]
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	# --- analyze data --- #
	output_file = output_folder + "shannon.txt"
	OG_sizes = []
	OG_shis = []
	OG_eveness = []
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			out.write( "OrthoGroupID\tSize\tShannonIndex\tEveness\n" )	#header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				values_per_species = []
				for each in parts[1:]:
					if "," in each:
						values_per_species.append( each.count(',')+1 )
					elif len( each.strip() ) > 1:
						values_per_species.append( 1 )
					else:
						values_per_species.append( 0 )
				shi, eveness = calculate_shannon_index( values_per_species )
				OG_sizes.append( math.log( sum( values_per_species ), 10 ) )
				OG_shis.append( shi )
				OG_eveness.append( eveness )
				out.write( "\t".join( map( str, [ 	parts[0],	#OG ID
																	sum( values_per_species ),	#number of sequences
																	shi,	#shannon index
																	eveness#EVENESS is based on all represented species. Would it be better to consider all species in total dataset?
																] ) ) + "\n" )
				line = f.readline()
	
	# --- construct figure --- #
	fig_file = output_folder + "shannon.pdf"
	construct_figure( fig_file, OG_sizes, OG_shis, "orthogroup diversity [Shannon Index]" )
	
	# --- construct figure --- #
	fig_file = output_folder + "eveness.pdf"
	construct_figure( fig_file, OG_sizes, OG_eveness, "orthogroup diversity [Eveness]" )
	
	


if '--out' in sys.argv and '--in' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

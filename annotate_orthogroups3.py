### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """ python annotate_orthogroups3.py
							--out <OUTPUT_FILE>
							--in <FREQUENCY_INPUT_FILE>
							
							optional:
							--anno <ANNOTATION_FILE>
							feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""

import os, sys, glob
from operator import itemgetter

# --- end of imports --- #


def load_anno( anno_file ):
	"""! @brief load annotation """
	
	anno = {}
	with open( anno_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno_string = ".".join( parts[1:-1] ).replace( parts[0]+".", "" )
			anno.update( { parts[0]: anno_string } )
			line = f.readline()
	return anno


def main( arguments ):
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if '--anno' in arguments:
		anno_file = arguments[ arguments.index('--anno')+1 ]
		anno = load_anno( anno_file )
	else:
		anno = {}
	
	n = 3
	
	# --- construct output file --- #
	with open( output_file, "w", 0 ) as out:
		out.write( "OrthoGroupID\tAnno\n" )
		with open( input_file, "r" ) as f:
			f.readline()	#remove header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if "; " in parts[1]:
					subparts = parts[1].split('; ')
					counts = []
					for hit in list(set( subparts )):
						counts.append( { 'hit': hit, 'count': subparts.count( hit ) } )
					best_hits = sorted( counts, key=itemgetter('count') )[::-1]
					fin_hit = []
					for i in range( n ):
						try:
							try:
								fin_hit.append( best_hits[ i ]['hit'] + "(" + str( round( 100.0 * best_hits[ i ]['count'] / len( subparts ), 2 ) ) + "%). " + anno[ best_hits[ i ]['hit'] ] )
							except KeyError:
								fin_hit.append( best_hits[ i ]['hit'] + "(" + str( round( 100.0 * best_hits[ i ]['count'] / len( subparts ), 2 ) ) + "%). n/a" )
						except IndexError:
							pass
					annotation = "; ".join( fin_hit )
				else:
					annotation = parts[1]
				out.write( "\t".join( [ parts[0], annotation ] ) + "\n" )
				line = f.readline()


if '--out' in sys.argv and '--in' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

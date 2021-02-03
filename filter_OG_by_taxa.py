### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """ python filter_OG_by_taxa.py
							--out <OUTPUT_OG_FILE>
							--in <INPUT_OG_FILE>
							--taxa <MINIMAL_TAXA_NUMBER> | --optimize
							feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
			"""

import re, os, sys
from operator import itemgetter


# --- end of imports --- #



def main( arguments ):
	""""! @brief run everything """
	
	output_file = arguments[ arguments.index('--out')+1 ]
	input_file = arguments[ arguments.index('--in')+1 ]
	if "--optimize" in arguments:
		optimization = True
		taxas = range( 100 )
	else:
		taxas = [ int( arguments[ arguments.index('--taxa')+1 ] ) ]
	
	for taxa in taxas:
		kept_OGs = 0
		dropped_OGs = 0
		with open( output_file, "w" ) as out:
			with open( input_file, "r" ) as f:
				out.write( f.readline() )	#header
				line = f.readline()
				while line:
					parts = line.strip().split('\t')
					counter = 0
					for each in parts[1:]:
						if len( each.strip() ) > 1:
							counter += 1
					if counter >=taxa:
						out.write( line )
						kept_OGs += 1
					else:
						dropped_OGs += 1
					line = f.readline()
		
		print "number of surviving OGs: " + str( kept_OGs ) + "\tnumber of dropped OGs: " + str( dropped_OGs ) + "\tcutoff: " + str( taxa )


if '--out' in sys.argv and '--in' in sys.argv and '--taxa' in sys.argv:
	main( sys.argv )
elif '--out' in sys.argv and '--in' in sys.argv and '--optimize' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

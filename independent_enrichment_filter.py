### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python independent_enrichment_filter.py
					--in <INPUT_FILES(COMMA-SEPARAED_LIST)>
					--out <OUTPUT_FILE>
					
					optional:
					--ecut <MIN_ENRICHMENT_CUTOFF>[1]
					--ncut <MIN_OCCUPANCY_CUTOFF>[all]
					--names <NAMES_OF_INPUT_SAMPLES>[sampleX]
					"""


import os, sys
from operator import itemgetter

# --- end of imports --- #

def load_og_enrichment_data( input_file, anno, ecutoff ):
	"""! @brief load orthogroup enrichment result data """
	
	data = []
	with open( input_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				anno[ parts[0] ]
			except KeyError:
				anno.update( { parts[0]: parts[5] } )
			if float( parts[2] ) > ecutoff:	#only consider enriched orthogroups
				data.append( { 'og': parts[0], 'enrichment': float( parts[2] ) } )
			line = f.readline()
	return data


def main( arguments ):
	"""! @brief run everything """
	
	input_files = arguments[ arguments.index('--in')+1 ].split(',')
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if '--names' in arguments:
		names = arguments[ arguments.index('--names')+1 ].split(',')
	else:
		names = []
		for i in range( len( input_files ) ):
			names.append( "sample" + str( i+1 ) )
	
	if '--ncut' in arguments:
		ncutoff = float( arguments[ arguments.index('--ncut')+1 ] )
	else:
		ncutoff = len( input_files )	#needs to be present in all datasets
	
	if '--ecut' in arguments:
		ecutoff = float( arguments[ arguments.index('--ecut')+1 ] )
	else:
		ecutoff = 1.0	#needs to be present in all datasets
	
	
	
	# --- load data --- #
	anno = {}
	collected_data = []
	for filename in input_files:
		collected_data.append( load_og_enrichment_data( filename, anno, ecutoff ) )
	
	# --- check for consistent enrichment in lineages --- #
	enrichment = {}
	for idx, data in enumerate( collected_data ):
		for og in data:
			try:
				enrichment[ og['og'] ]['vals'].append( og['enrichment'] )
				enrichment[ og['og'] ]['names'].append( names[ idx ] )
			except KeyError:
				enrichment.update( { og['og']: { 'vals': [ og['enrichment'] ], 'names': [ names[ idx ] ] } } )
	
	fin_data = []
	occupancy_distribution = []
	for og in enrichment.keys():
		occupancy_distribution.append( len( enrichment[ og ]['vals'] ) )
		if len( enrichment[ og ]['vals'] ) >= ncutoff:
			fin_data.append( { 	'og': og,
												'avg': sum( enrichment[ og ]['vals'] ) / len(enrichment[ og ]['vals']  ),
												'vals': ";".join( map( str, enrichment[ og ]['vals'] ) ),
												'names': ";".join( enrichment[ og ]['names'] ),
												'anno': anno[ og ]
										} )
	
	# --- generate output file --- #
	with open( output_file, "w" ) as out:
		out.write( "OrthoGroupID\tAverageEnrichment\tValues\tLineages\tAnnotation\n" )
		for entry in sorted( fin_data, key=itemgetter('avg') )[::-1]:
			out.write( "".join( map( str, [ entry['og'], entry['avg'], entry['vals'], entry['names'], entry['anno'] ] ) ) + "\n" )
	
	# --- statitstics --- #
	for i in range( len( input_files ) ):
		print "number of orthogroups with enrichment in " + str( i+1 ) + " lineage(s): " + str( occupancy_distribution.count( i+1 ) )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

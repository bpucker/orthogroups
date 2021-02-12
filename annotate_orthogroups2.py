### Boas Pucker ###
### v0.1 ###
### bpucker@cebitec.uni-bielefeld.de ###


__usage__ = """ python annotate_orthogroups2.py
							--out <OUTPUT_FILE>
							--og <ORTHOGROUP_FILE>
							--in <BLAST_RESULT_FILE>
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
			anno_string = ".".join( parts[:-1] ).replace( parts[0]+"."+parts[0]+".", parts[0]+"." )
			anno.update( { parts[0]: anno_string } )
			line = f.readline()
	return anno


def load_BLAST_results( blast_result_file, n=3 ):
	"""! @brief load top3 annotations per orthogroup """
	
	# --- get top hit per query sequence --- #
	best_hit = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if float( parts[-1] ) > best_hit[ parts[0] ]:
					best_hit[ parts[0] ] = { 'score': float( parts[-1] ), 'hit': parts[1] }
			except KeyError:
				best_hit.update( { parts[0]: { 'score': float( parts[-1] ), 'hit': parts[1] } } )
			line = f.readline()
	
	print "number of initial BLAST hits: " + str( len( best_hit.keys() ) )
	
	# --- combine hits per orthogroup --- #
	best_hits_per_og = {}
	for key in best_hit.keys():
		try:
			best_hits_per_og[ key.split('_%_')[0] ].append( best_hit[key]['hit'] )
		except KeyError:
			best_hits_per_og.update( { key.split('_%_')[0]: [ best_hit[key]['hit'] ] } )
	
	print "number of orthogroups with BLAST hits: " + str( len( best_hits_per_og.keys() ) )
	
	return best_hits_per_og
	# # --- get the n most frequent hits --- #
	# final_hits_per_og = {}
	# for key in best_hits_per_og.keys():
		# hits = best_hits_per_og.keys()
		# counts = []
		# for hit in list(set( hits )):
			# counts.append( { 'hit': hit, 'count': hits.count( hit ) } )
		# best_hits = sorted( counts, key=itemgetter('count') )[::-1]
		# fin_hit = []
		# for i in range( n ):
			# try:
				# fin_hit.append( best_hits[ i ]['hit'] )
			# except IndexError:
				# pass
		# final_hits_per_og.update( { key: fin_hit } )
	# return final_hits_per_og


def main( arguments ):
	
	blast_result_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	orthogroup_file = arguments[ arguments.index('--og')+1 ]
	
	if '--anno' in arguments:
		anno_file = arguments[ arguments.index('--anno')+1 ]
		anno = load_anno( anno_file )
	else:
		anno = {}
	
	
	# --- load BLAST results --- #
	blast_results = load_BLAST_results( blast_result_file )
	print "number of annotations loaded for orthogroups: " + str( len( blast_results.keys() ) )
	
	# --- construct output file --- #
	with open( output_file, "w" ) as out:
		out.write( "OrthoGroupID\tAnno\n" )
		with open( orthogroup_file, "r" ) as f:
			f.readline()	#remove header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				new_line = [ parts[0] ]
				try:
					annotation = []
					for hit in sorted( list( set( blast_results[ parts[0] ] ) ) ):
						try:
							annotation.append( anno[ hit ] )
						except KeyError:
							annotation.append( hit )
					new_line.append( "; ".join( annotation ) )
				except KeyError:
					pass
				while len( new_line ) < 2:
					new_line.append( "." )
				out.write( "\t".join( new_line ) + "\n" )
				line = f.readline()
	
	# --- construct frequency output file --- #
	with open( output_file+".freq", "w" ) as out:
		out.write( "OrthoGroupID\tAnnoFreq\n" )
		with open( orthogroup_file, "r" ) as f:
			f.readline()	#remove header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				new_line = [ parts[0] ]
				try:
					new_line.append( "; ".join( blast_results[ parts[0] ] ) )
				except KeyError:
					pass
				while len( new_line ) < 2:
					new_line.append( "." )
				out.write( "\t".join( new_line ) + "\n" )
				line = f.readline()

if '--out' in sys.argv and '--in' in sys.argv and '--og' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

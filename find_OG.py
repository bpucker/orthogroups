### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python find_OG.py
					--in <SEQ_INPUT_FILE>
					--ref <OG_SEQ_FASTA>
					--out <OUTPUT_FOLDER>
					"""


import os, sys
from operator import itemgetter

# --- end of imports --- #

def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if float( parts[-1] ) > best_hits[ parts[0] ]['score']:
					best_hits[ parts[0] ] = { 'score': float( parts[-1] ), 'hit': parts[1] }
			except:
				best_hits.update( { parts[0]: { 'score': float( parts[-1] ), 'hit': parts[1] } } )
			line = f.readline()
	return best_hits



def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	ref_file = arguments[ arguments.index('--ref')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	cpus = 4
	
	if output_folder[ -1 ] != "/":
		output_folder += "/"
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	blast_db = output_folder + "blastdb"
	os.popen( "makeblastdb -in " + ref_file + " -out " + blast_db + " -dbtype prot" )
	
	blast_result_file = output_folder + "blast_results.txt"
	if not os.path.isfile( blast_result_file ):
		os.popen( "blastp -query " + input_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.00001 -num_threads " + str( cpus ) )
	
	# --- get best BLAST match per query --- #
	best_match = load_best_blast_hit( blast_result_file )
	best_hit_file = output_folder + "best_hit.txt"
	with open( best_hit_file, "w" ) as out:
		out.write( "Query\tSubject\tScore\n" )
		for key in best_match:
			out.write( "\t".join( map( str, [ key, best_match[key]['hit'], best_match[key]['score'] ] ) ) + "\n" )
	
	# --- find best matching OGs --- #
	OGs = []
	for hit in best_match.values():
		OGs.append( hit['hit'].split('_%_')[0] )
	
	output_data = []
	for OG in list( set( OGs ) ):
		output_data.append( { 'og': OG, 'counts': OGs.count( OG ) } )
	
	fin_doc_file = output_folder + "best_OG_hit.txt"
	with open( fin_doc_file, "w" ) as out:
		out.write( "OG\tNumberOfHits\n" )
		for entry in sorted( output_data, key=itemgetter('counts') )[::-1]:
			out.write( entry['og'] + "\t" + str( entry['counts'] ) + "\n" )


if '--in' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

require 'rubygems'
require 'bio'
 
# creating a Bio::Sequence::NA object containing ambiguous alphabets
#as = Bio::Sequence::NA.new("atgcyrwskmbdhvn")

#print as.to_s

#print "\n"


#seq1 = Bio::Sequence::AA.new("gggggg")
#seq2 = Bio::Sequence::AA.new("ggggt")
#seq3 = Bio::Sequence::AA.new("ggt")



seq1 = "KMLFGVVFFFGG"
seq2 ="LMGGHHF"
seq3 = "GKKKKGHHHGHRRRGR"
seq4 = "KKKGHHHGHRERGR"
seqs = [ seq1, seq2, seq3, seq4 ] 


# MAFFT
options = [ '--maxiterate', '1000', '--localpair' ]
mafft = Bio::MAFFT.new('/home/zma/SOFTWARE/mafft-6.847-without-extensions/scripts/mafft', options )
report = mafft.query_align( seqs)
#puts report.alignment.output_fasta
report.alignment.each { |x| puts x.to_s }
puts 'OK'
puts

#clustalw
clustalw = Bio::ClustalW.new('/home/zma/SOFTWARE/clustalw-2.1/src/clustalw2' )
report = clustalw.query_align( seqs)
#puts report.alignment.output_fasta.to_s
report.alignment.each { |x| puts x.to_s }
puts 'OK'
puts

#muscle
options = [ '-quiet', '-maxiters', '64' ]
muscle = Bio::Muscle.new('/home/zma/SOFTWARE/muscle3.8.31/src/muscle', options )
report = muscle.query_align( seqs)
#puts report.alignment.output_fasta.to_s
report.alignment.each { |x| puts x.to_s }
puts 'OK'
puts
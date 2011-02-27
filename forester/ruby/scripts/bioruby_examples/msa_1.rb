require 'rubygems'
require 'bio'
 
#############

# Reads in a clustalw formatted multiple sequence alignment
# from a file named "infile_clustalw.aln" and stores it in 'report'.
report = Bio::ClustalW::Report.new(File.read('infile_clustalw.aln'))

# Accesses the actual alignment.
align = report.alignment

# Goes through all sequences in 'align' and prints the
# actual molecular sequence.
align.each do |entry|
  puts entry.seq
end


# Creates a new file named "outfile.fasta" and writes
# multiple sequence alignment 'align' to it in fasta format.
File.open('outfile.fasta', 'w') do |f|
  f.write(align.output(:fasta))
end

# Creates a new file named "outfile.aln" and writes
# multiple sequence alignment 'align' to it in clustal format.
File.open('outfile.aln', 'w') do |f|
  f.write(align.output(:clustal))
end

#############

seq1 = Bio::Sequence.auto("gggggg")


puts seq1.output(:fasta)
#seq2 = Bio::Sequence::AA.new("ggggt")
#seq3 = Bio::Sequence::AA.new("ggt")



seqs = [ "MFQIPEFEPSEQEDSSSAER",
         "MGTPKQPSLAPAHALGLRKS",
         "PKQPSLAPAHALGLRKS",
         "MCSTSGCDLE" ] 


# MAFFT
options = [ '--maxiterate', '1000', '--localpair' ]
mafft = Bio::MAFFT.new('/home/zma/SOFTWARE/mafft-6.847-without-extensions/scripts/mafft', options )
report = mafft.query_align( seqs)

# Accesses the actual alignment
align = report.alignment

# Prints each sequence to the console.
align.each { |s| puts s.to_s }


puts 'MAFFT OK'
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

file = Bio::FastaFormat.open('bcl2.fasta')
file.each do |entry|
  puts entry.entry_id           # Gets the identifier, e.g. 'sp|O35147|BAD_RAT'.
  puts entry.definition         # Gets the complete fasta description line.
  puts entry.seq                # Gets the actual sequence.
  puts entry.aaseq.composition  # Gets the amino acid composition. 
end
puts 'OK'
puts

Bio::FlatFile.auto('bcl2.fasta') do |ff|
  ff.each do |entry|
    puts entry.entry_id           # Gets the identifier, e.g. 'sp|O35147|BAD_RAT'.
    puts entry.definition         # Gets the complete fasta description line.
    puts entry.seq                # Gets the actual sequence.
    puts entry.aaseq.composition  # Gets the amino acid composition.
  end
end
puts 'OK'
puts







require 'bio'

seq_ary = Array.new
ff = Bio::FlatFile.auto('bcl2.fasta')
ff.each_entry do |entry|
  seq_ary.push(entry)
  puts entry.entry_id          # prints the identifier of the entry
  puts entry.definition        # prints the definition of the entry
  puts entry.seq               # prints the sequence data of the entry
end

# Creates a multiple sequence alignment (possibly unaligned) named
# 'seqs' from array 'seq_ary'.
seqs = Bio::Alignment.new( seq_ary )


seqs.each { |seq| puts seq.to_s }


puts seqs.consensus

# Writes multiple sequence alignment (possibly unaligned) 'seqs'
# to a file in PHYLIP format.
File.open('out0.phylip', 'w') do |f|
  f.write(seqs.output(:phylip))
end

File.open('out0.fasta', 'w') do |f|
  f.write(seqs.output(:fasta))
end

exit
#############

# Reads in a FASTA-formatted multiple sequence alignment (which does
# not have to be aligned, though) and stores its sequences in
# array 'seq_ary'.
seq_ary = Array.new
fasta_seqs = Bio::Alignment::MultiFastaFormat.new(File.open('bcl2.fasta').read)
fasta_seqs.entries.each do |seq|
  seq_ary.push( seq )
end

# Creates a multiple sequence alignment (possibly unaligned) named
# 'seqs' from array 'seq_ary'.
seqs = Bio::Alignment.new( seq_ary )
seqs.each { |seq| puts seq.to_s }


puts seqs.consensus

# Writes multiple sequence alignment (possibly unaligned) 'seqs'
# to a file in PHYLIP format.
File.open('out1.phylip', 'w') do |f|
  f.write(seqs.output(:phylip))
end

File.open('out1.fasta', 'w') do |f|
  f.write(seqs.output(:fasta))
end

exit
#################

#ff = Bio::FlatFile.new(Bio::FastaFormat, 'bcl2.fasta')
#ff.each_entry do |f|
#  puts "definition : " + f.definition
#  puts "nalen      : " + f.nalen.to_s
#  puts "naseq      : " + f.naseq
#end
#exit

seq_ary = Array.new
Bio::FastaFormat.open('bcl2.fasta') do | file |
  file.each do |entry|
    puts entry.entry_id           # Gets the identifier, e.g. 'sp|O35147|BAD_RAT'.
    # puts entry.definition         # Gets the complete fasta description line.
    #puts entry.seq                # Gets the actual sequence.
    seq_ary.push( entry )
  end
end
seqs =Bio::Alignment.new( seq_ary )
seqs.each { |x| puts x }
puts seqs.consensus
exit



#will be obsolete!
#seqs =Bio::Alignment.readfiles(File.open('bcl2.fasta'))
#seqs.entries.each do |seq|
#  puts seq.to_biosequence
#end

#Bio::Alignment.



#exit
#############

seqs = Bio::Alignment::MultiFastaFormat.new(File.open('bcl2.fasta').read)
seqs.entries.each do |seq|
  puts seq.to_seq.output(:genbank)
end
puts
puts
puts :genbank
puts seqs.entries[0].to_seq.output(:genbank)
puts
puts :fasta
puts seqs.entries[0].to_seq.output(:fasta)
puts
puts :embl
puts seqs.entries[0].to_seq.output(:embl)
puts
puts :raw
puts seqs.entries[0].to_seq.output(:raw)
puts
puts :fasta_ncbi
puts seqs.entries[0].to_seq.output(:fasta_ncbi)
puts
puts :fastq
puts seqs.entries[0].to_seq.output(:fastq)
puts
puts :fastq_sanger
puts seqs.entries[0].to_seq.output(:fastq_sanger)
puts
puts :fastq_solexa
puts seqs.entries[0].to_seq.output(:fastq_solexa)
puts
puts :fastq_illumina
puts seqs.entries[0].to_seq.output(:fastq_illumina)
puts
puts :fasta_numeric
puts seqs.entries[0].to_seq.output(:fasta_numeric)
puts
puts :qual
puts seqs.entries[0].to_seq.output(:qual)
#exit
##############


# Reads in a ClustalW formatted multiple sequence alignment
# from a file named "infile_clustalw.aln" and stores it in 'report'.
report = Bio::ClustalW::Report.new(File.read('infile_clustalw.aln'))

# Accesses the actual alignment.
msa = report.alignment

# Goes through all sequences in 'msa' and prints the
# actual molecular sequence.
msa.each do |entry|
 # puts entry.seq
end

##############

DEFAULT_PARSER = Bio::Alignment::MultiFastaFormat
puts DEFAULT_PARSER.to_s

#file = Bio::Alignment.readfiles('bcl2.fasta', Bio::Alignment::MultiFastaFormat)
#file.each do |entry|
#  puts entry.entry_id           # Gets the identifier, e.g. 'sp|O35147|BAD_RAT'.
#  puts entry.definition         # Gets the complete fasta description line.
#  puts entry.seq                # Gets the actual sequence.
  #puts entry.aaseq.composition  # Gets the amino acid composition.
#end
#puts 'OK'
#puts

file = Bio::FastaFormat.open('bcl2.fasta')
file.each do |entry|
   puts entry.entry_id           # Gets the identifier, e.g. 'sp|O35147|BAD_RAT'.
  puts entry.definition         # Gets the complete fasta description line.
  puts entry.seq                # Gets the actual sequence.
  # do something on each fasta sequence entry
end


##############

# Creates a new file named "outfile.fasta" and writes
# multiple sequence alignment 'msa' to it in fasta format.
File.open('outfile.fasta', 'w') do |f|
  f.write(msa.output(:fasta))
end

# Other formats
File.open('outfile.clustal', 'w') do |f|
  f.write(msa.output(:clustal))
end
File.open('outfile.phylip', 'w') do |f|
  f.write(msa.output(:phylip))
end
File.open('outfile.phylipnon', 'w') do |f|
  f.write(msa.output(:phylipnon))
end
File.open('outfile.msf', 'w') do |f|
  f.write(msa.output(:msf))
end
File.open('outfile.molphy', 'w') do |f|
  f.write(msa.output(:molphy))
end



#############

seq1 = Bio::Sequence.auto("gggggg")


puts seq1.output(:fasta)
#seq2 = Bio::Sequence::AA.new("ggggt")
#seq3 = Bio::Sequence::AA.new("ggt")



seqs = ['MFQIPEFEPSEQEDSSSAER',
        'MGTPKQPSLAPAHALGLRKS',
        'PKQPSLAPAHALGLRKS',
        'MCSTSGCDLE'] 


# MAFFT
options = [ '--maxiterate', '1000', '--localpair' ]
mafft = Bio::MAFFT.new('mafft', options )
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







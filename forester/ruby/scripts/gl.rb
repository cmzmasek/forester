File.open( "gl.txt" ) do | file |
  while line = file.gets
    if line =~ /^([A-Z0-9]{3,5})/
      puts $1 + "\t" + "/home/czmasek/DATA/SEQUENCES/GENOMES/PROTEIN_PREDICTIONS_CDHIT_095/ALL/" + $1 + ".fasta"
    else
      puts line
    end
  end
end
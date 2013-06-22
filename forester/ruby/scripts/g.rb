File.open( "g" ) do | file |
  while line = file.gets
    if line =~ /^(\S{3,5})\s+\S{3,5}/
      puts "/home/czmasek/DATA/GENOME_HMMSCAN/CDHIT_095_PFAM270X/" + $1 + ".hmmscan    " + $1
    else
      puts line
    end
  end
end
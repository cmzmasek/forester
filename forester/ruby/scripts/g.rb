File.open( "g" ) do | file |
  while line = file.gets
    if line =~ /^([A-Z0-9]{3,5})/
      puts "/home/czmasek/DATA/GENOME_HMMSCAN/CDHIT_095_PFAM270X/" + $1 + ".hmmscan    " + $1
    else
      puts line
    end
  end
end
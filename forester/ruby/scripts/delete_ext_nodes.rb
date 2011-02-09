#!/usr/local/bin/ruby -w

infile = ARGV[ 0 ]

metazoa_choanoflagellata = [ 
"Metazoa_Choanoflagellata", 
"Metazoa",
"Bilateria_Cnidaria",
"Bilateria",
"Deuterostomia",
"Chordata",
"Urochordata_Vertebrata",
"Vertebrata",
"Tetrapoda",
"Amniota",
"Eutheria",
"Euarchontoglires",
"Primates",
"Rodentia",
"Teleostei",
"Euteleostei",
"Smegmamorpha",
"Tetraodontiformes",
"Urochordata",
"Ascidiacea",
"Urochordata",
"Protostomia",
"Ecdysozoa",
"Arthropoda",
"Insecta",
"Lepidoptera_Diptera_Hymenoptera",
"Diptera",
"Culicoidea",
"Hymenoptera",
"Nematoda",
"Annelida_Mollusca",
"Annelida" ]

if infile == nil
 puts "no infile"
 exit
end

File.open( infile ) do | file |
    while line = file.gets
        if line =~ /^[0-9A-Z]{3,5}\s/
        elsif line =~ /^\t/
        elsif line =~ /^{/
         elsif line =~ /^f_\d/
        else   
            line =~ /(\S+)/
            first = $1
            if metazoa_choanoflagellata.include?( first ) 
                puts( line )
            end
        end
    end
end 

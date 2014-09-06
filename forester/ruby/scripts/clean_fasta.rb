require 'fileutils'

infile = ARGV[ 0 ]

File.open(infile, 'r') do |f|
   while line = f.gets
       if line =~ /(>.+)~.+/
          line = $1
       elsif line =~ /(>.+)\|NAME=.+/
           line = $1
       end
       puts line
   end
end




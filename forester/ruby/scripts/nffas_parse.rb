#!/usr/local/bin/ruby -w

infile = ARGV[ 0 ]

File.open( infile ) do | file |
  while line = file.gets
    if line =~ /^>>(.+)/
      puts ">" + $1
    elsif line =~ /^([A-Z]+)/
      puts $1
    end
  end
end

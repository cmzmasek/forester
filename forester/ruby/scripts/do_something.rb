require 'set'
require 'fileutils'

set = Set.new

File.open('infile', 'r') do |f|
   while line = f.gets
     set << line
   end
end

counter = 0
set.each do |s| 
  puts s
  s = s.chomp + '.prot'
  if File.file?( s ) && File.size?( s )
    FileUtils.cp( s, "keep/" + s )
    counter += 1 
  end
end

puts "Copied " + counter.to_s + " files."


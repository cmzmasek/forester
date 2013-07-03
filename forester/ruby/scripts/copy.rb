require 'fileutils'
Dir.foreach(".") do |f|
  next if f == '.' or f == '..'
  if f =~ /^([A-Z0-9]{3,5}).*\..+/
    tax = $1
    FileUtils.cp f, "all/" + tax + ".fasta"
  elsif 
    puts "ignored " + f.to_s
  end
end
PRG = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/hsp.rb"
Dir.foreach(".") do |f|
  next if f == '.' or f == '..'
  if f =~ /^([A-Z0-9]{3,6})\.hmmscan_260/
    taxon = $1
    cmd = "#{PRG} -ie=10 -a -pe=1 -m=RRMa -s=#{taxon} #{taxon}.hmmscan_260 #{taxon}_ie10_domain_table > #{taxon}_ie10_summary_table.txt"
    puts cmd
    output = IO.popen( cmd )
  elsif f =~ /^([A-Z0-9]{3,6})_cdhit090\.hmmscan_260/
    taxon = $1
    cmd = "#{PRG} -ie=10 -a -pe=1 -m=RRMa -s=#{taxon} #{taxon}_cdhit090.hmmscan_260 #{taxon}_cdhit090_ie10_domain_table > #{taxon}_cdhit090_ie10_summary_table.txt"
    puts cmd
    output = IO.popen( cmd )
  end
end
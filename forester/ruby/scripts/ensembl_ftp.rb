require 'net/ftp'

EMAIL           = '?'
PUB_RELEASE_DIR = '/pub/release-65/fasta'
PEP_DIR         = '/pep'

ftp = Net::FTP.new('ftp.ensembl.org', 'anonymous', EMAIL)
ftp.passive = true # To avoid "No route to host" error.
ftp.chdir( PUB_RELEASE_DIR )
files = ftp.list('*_*') # To only look at files with an underscore.
count = 0
files.each do | file |
  species = file.split().last
  begin
    ftp.chdir(species + PEP_DIR)
    pepfiles = ftp.list()
    pepfiles.each do | pepfile |
      pepfile = pepfile.split().last
      if pepfile =~ /all.fa.gz/ # Only want the "all.fa.gz" files (and not the
                                # "abinitio" files).
        ftp.getbinaryfile(pepfile)
        puts 'downloaded "' + pepfile + '"'
        count += 1
      end
    end
  rescue Exception
    puts 'ignoring "' + species + '"'
  end
  ftp.chdir(PUB_RELEASE_DIR) # To go back to the starting directory.
end
ftp.close
puts 'done (downloaded ' + count.to_s + ' files)'
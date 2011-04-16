
NEW_DIR = '/home/czmasek/WORK/GENOME_HMMPFAM/HMMSCAN3_no_bias_250/'
OLD_DIR = '/home/czmasek/WORK/GENOME_HMMPFAM/HMMSCAN30b3_no_bias_240/'

Dir.foreach('.') { |f| 
  if File.symlink?( f )
     link = File.readlink( f )
     puts f + ' -> ' + link 
   
     link =~ /\/([^\/]+)\/([^\/]+)\.hmmscan30b3_240/
     group = $1
     core_name = $2
     puts '  => ' + group + ' / ' + core_name
     new_link = NEW_DIR + group + '/' + core_name + '.hmmscan_250'
     puts '  => ' + new_link
     puts
     if ( !File.exists?( new_link ) || !File.exists?( link ) )
        puts 'ERROR!'
        exit
     end
     #what is called new_link will by linked to by 'newdir/' + f.to_s:
     File.symlink( new_link, 'newdir/' + f.to_s )
  end
}


#File.symlink("testfile", "link2test") 

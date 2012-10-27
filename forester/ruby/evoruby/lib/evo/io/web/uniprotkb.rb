
require 'net/http'
require 'uri'

require 'lib/evo/io/parser/uniprot_parser'

module Evoruby


  class UniprotKB

    BASE_URL = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb;"


    def UniprotKB::get_by_id( id, style = "raw", format = "uniprot" )
      url_str = BASE_URL + "id=#{id};format=#{format};style=#{style}"
      uri = URI.parse url_str
      response = Net::HTTP.get_response uri
      lines = []
      response.body.each_line do |line|
        lines << line
        puts line
      end
      lines
    end


    def UniprotKB::get_entry_by_id( id  )
      lines = get_by_id( id,  "raw", "uniprot" )
      p = UniprotParser.new
      return p.parse( lines )
    end

  end

end

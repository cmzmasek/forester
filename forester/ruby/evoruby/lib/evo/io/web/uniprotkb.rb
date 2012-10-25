
require 'net/http'
require 'uri'

module Evoruby


  class UniprotKB
    def initialize

    end

    def get
      require 'net/http'
      require 'uri'

      uri = URI.parse("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb;id=1433X_MAIZE;format=uniprot;style=raw")
      response = Net::HTTP.get_response uri
      puts response.body

    end

  end

end

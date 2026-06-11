# forester -- software libraries and applications
# for evolutionary biology and genomics.
# Copyright (C) 2026 Christian M. Zmasek
# All rights reserved
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact: czmasek at jcvi dot org

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

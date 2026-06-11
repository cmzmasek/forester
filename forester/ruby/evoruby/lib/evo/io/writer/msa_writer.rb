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

require 'lib/evo/util/constants'
require 'lib/evo/util/util'

module Evoruby

  class MsaWriter

    def initialize()
      raise TypeError, "Cannot instanciate abstract class MsaWriter"
    end

    def set_max_name_length( length )
    end

    def set_exception_if_name_too_long( exception_if_name_too_long )
    end

    def write( msa, path )
    end

  end # class MsaWriter

end # module Evoruby

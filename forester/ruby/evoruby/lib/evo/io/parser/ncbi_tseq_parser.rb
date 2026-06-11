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

require 'lib/evo/io/parser/msa_parser'
require 'lib/evo/taxonomy/taxonomy'
require 'lib/evo/msa/msa'

require 'iconv'

module Evoruby

    class NcbiTSeqParser < MsaParser

        TSEQ_SEQ = "TSeq_sequence"
        TSEQ_DEFLINE = "TSeq_defline"
        TSEQ_ORGNAME = "TSeq_orgname"
        TSEQ_TAXID = "TSeq_taxid"
        TSEQ_SID = "TSeq_sid"
        TSEQ_ACCVER = "TSeq_accver"
        TSEQ_GI = "TSeq_gi"
        TSEQ_TYPE = "TSeq_seqtype"
        TSEQ_LENGTH = "TSeq_length"

        def initialize
        end


        #  <TSeqSet>
        #<TSeq>
        #  <TSeq_seqtype value="protein"/>
        #  <TSeq_gi>29341016</TSeq_gi>
        #  <TSeq_accver>AAO78806.1</TSeq_accver>
        #  <TSeq_sid>gnl|mbpwusl|BT3701</TSeq_sid>
        #  <TSeq_taxid>226186</TSeq_taxid>
        #  <TSeq_orgname>Bacteroides thetaiotaomicron VPI-5482</TSeq_orgname>
        #  <TSeq_defline>SusD [Bacteroides thetaiotaomicron VPI-5482]</TSeq_defline>
        #  <TSeq_length>551</TSeq_length>
        #  <TSeq_sequence>MKTKYIKQLFSAALIAVLSSGVTSCINDLDISPIDPQTGGSFDQQGVFVKGYAMLGVTGQKGIDGSPDLDGQDEGESGFYRTTFNCNELPTDECLWAWQENQDIPQLTSISWSPSSQRTEWVYVRLGYDITQYNFFLDQTEGMTDAETLRQRAEIRFLRALHYWYFLDLFGKAPFKEHFSNDLPVEKKGTELYTYIQNELNEIEADMYEPRQAPFGRADKAANWLLRARLYLNAGVYTGQTDYAKAEEYASKVIGSAYKLCTNYSELFMADNDENENAMQEIILPIRQDGVKTRNYGGSTYLVCGTRVAGMPRMGTTNGWSCIFARAAMVQKFFSNLEDVPMLPADVEIPTKGLDTDEQIDAFDAEHGIRTEDMIKAAGDDRALLYSGVGGGRRKIQTDAISGFTDGLSIVKWQNYRSDGKPVSHATYPDTDIPLFRLAEAYLTRAEAIFRQGGDATGDINELRKRANCTRKVQTVTEQELIDEWAREFYLEGRRRSDLVRFGMFTTNKYLWDWKGGAMNGTSVASYYNKYPIPVSDINNNRNMSQNEGYK</TSeq_sequence>
        #</TSeq>

        def parse( path )
            Util.check_file_for_readability( path )
            seqs = Msa.new

            in_seq        = false
            gi = nil
            accver = nil
            sid = nil
            taxid = nil
            orgname = nil
            defline = nil
            seq_str = nil
            line_counter = 1
            ic = Iconv.new( 'UTF-8//IGNORE', 'UTF-8' )
            File.open( path ) do | file |
                while line = file.gets
                    line = ic.iconv( line )
                    line_counter += 1
                    if can_ignore?( line )

                    elsif line =~ /^\s*<TSeq>/
                        in_seq = true


                    elsif in_seq
                        if line =~ /^\s*<\/TSeq>/
                            in_seq = false
                            taxonomy = nil
                            if taxid != nil || orgname != nil
                                id_source = nil
                                if taxid != nil
                                    id_source = "ncbi"
                                end
                                taxonomy = Taxonomy.new( orgname, taxid , id_source )
                            end
                            id = nil
                            id_source = nil
                            symbol = nil
                            if gi != nil
                                id = gi
                                id_source = "gi"
                                if sid != nil
                                    symbol = sid
                                elsif accver != nil
                                    symbol = accver
                                end
                            elsif sid != nil
                                id = sid
                                id_source = "ncbi"
                                if accver != nil
                                    symbol = accver
                                end
                            elsif accver != nil
                                id = accver
                                id_source = "ncbi"
                            end

                            sequence = Sequence.new( defline,
                                seq_str,
                                id,
                                id_source,
                                taxonomy,
                                symbol )

                            seqs.add_sequence( sequence )
                            gi = nil
                            accver = nil
                            sid = nil
                            taxid = nil
                            orgname = nil
                            defline = nil
                            seq_str = nil
                        elsif line =~ /^\s*<#{TSEQ_GI}>(\d+)<\/#{TSEQ_GI}>/
                            gi = $1
                        elsif line =~ /^\s*<#{TSEQ_ACCVER}>(.+)<\/#{TSEQ_ACCVER}>/
                            accver = $1
                        elsif line =~ /^\s*<#{TSEQ_SID}>(.+)<\/#{TSEQ_SID}>/
                            sid = $1
                        elsif line =~ /^\s*<#{TSEQ_TAXID}>(\d+)<\/#{TSEQ_TAXID}>/
                            taxid = $1
                        elsif line =~ /^\s*<#{TSEQ_ORGNAME}>(.+)<\/#{TSEQ_ORGNAME}>/
                            orgname = $1
                        elsif line =~ /^\s*<#{TSEQ_DEFLINE}>(.+)<\/#{TSEQ_DEFLINE}>/
                            defline = $1
                        elsif line =~ /^\s*<#{TSEQ_SEQ}>(.+)<\/#{TSEQ_SEQ}>/
                            seq_str = $1
                        elsif line =~ /^\s*<#{TSEQ_TYPE}/
                        elsif line =~ /^\s*<#{TSEQ_LENGTH}/
                        else
                            error_msg = "unexpected line format at line #{line_counter}: " + line
                            raise IOError, error_msg
                        end
                    end
                end
            end
            return seqs
        end

        private

        def can_ignore?( line )
            return ( line !~ /\S/ )
        end

    end # class NcbiTSeqParser

end # module Evoruby
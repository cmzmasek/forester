#!/usr/local/bin/ruby -w
#
# = exe/test - Test class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: test.rb,v 1.18 2010/10/08 22:04:17 cmzmasek Exp $
#
# last modified: 05/15/2007


require 'lib/evo/util/constants'
require 'lib/evo/taxonomy/taxonomy'
require 'lib/evo/sequence/sequence'
require 'lib/evo/msa/msa'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/sequence/domain_structure'
require 'lib/evo/sequence/protein_domain'
require 'lib/evo/table/basic_table'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/writer/phylip_sequential_writer'
require 'lib/evo/io/writer/nexus_writer'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/parser/ncbi_tseq_parser'
require 'lib/evo/io/parser/hmmsearch_domain_extractor'
require 'lib/evo/tool/domain_sequence_extractor'
require 'lib/evo/tool/hmmscan_summary'
require 'lib/evo/tool/domains_to_forester'
require 'lib/evo/io/parser/general_msa_parser'
require 'lib/evo/io/parser/basic_table_parser'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/soft/fastme'
require 'lib/evo/soft/tree_puzzle'



module Evoruby

    class Test

        GENERAL_MSA_FILE = "files/test/general_msa_file.txt"
        FASTA_FILE       = "files/test/fasta_file.txt"
        TSEQ_FILE        = "files/test/ncbi_tseq.xml"

        def initialize()
            @failures  = 0
            @successes = 0
        end



        def test_taxonomy()
            begin
                tax = Taxonomy.new( "pig" )

                if tax.get_name != "pig"
                    return false
                end

                tax1 = Taxonomy.new( "dog", "id", "source" )
                tax2 = tax1.copy

                if tax2.get_name != "dog"
                    return false
                end
                if tax2.get_id != "id"
                    return false
                end
                if tax2.get_id_source != "source"
                    return false
                end

                if !( tax1 == tax2 )
                    return false
                end

                if !( tax1 == tax1 )
                    return false
                end

                tax3 = Taxonomy.new( "dog", "id"  )
                if ( tax1 == tax3 )
                    return false
                end

                tax4 = Taxonomy.new( "dog" )
                tax5 = Taxonomy.new( "dog" )
                if !( tax4 == tax5 )
                    return false
                end

            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end


        def test_sequence()
            begin
                seq = Sequence.new( "seq1", "WLIQ" )
                if ( seq.get_length != 4 )
                    return false
                end
                if ( seq.get_residue( 3 ) != "Q" )
                    return false
                end
                seq.append!( "E?-*X_Y" )
                if ( seq.get_length != 11 )
                    return false
                end
                if ( seq.get_residue( 3 ) != "Q" )
                    return false
                end
                if ( seq.get_residue( 4 ) != "E" )
                    return false
                end
                seq.append!( "A V_" )
                if ( seq.get_length != 15 )
                    return false
                end
                if ( !Test::same?( seq.get_gap_length, 5 ) )
                    return false
                end
                if ( !Test::same?( seq.get_gap_ratio, 5.0 / 15.0 ) )
                    return false
                end
                seq.delete_residue!( 0 )
                seq.delete_residue!( 2 )
                seq2 = seq.copy()
                seq.delete_residue!( 0 )
                seq.delete_residue!( 0 )
                seq = nil
                if ( seq2.get_length != 13 )
                    return false
                end
                if ( seq2.get_sequence_as_string != "LIE?-*X_YA V_" )
                    return false
                end
                if ( seq2.get_slice( 2, 2 ) != "E?" )
                    return false
                end
                if ( seq2.get_slice( 0, 1 ) != "L" )
                    return false
                end
                if ( seq2.get_subsequence( 1, 4 ).get_sequence_as_string != "IE?-" )
                    return false
                end
                if ( seq2.get_name() != "seq1" )
                    return false
                end
                if ( seq2.get_slice!( 2, 2 ) != "E?" )
                    return false
                end
                if ( seq2.get_sequence_as_string != "LI-*X_YA V_" )
                    return false
                end
                if ( seq2.get_length != 11 )
                    return false
                end
                if ( seq2.get_character_code( 0 ) != 76 )
                    return false
                end
                str_0 = " Li-*X_YA V_ 3 3    1212 ?? B1J OU.Z "
                if ( Util.clean_seq_str( str_0 ) != "LI-X-YAV-XXXXXX-X" )
                    return false
                end

                tax = Taxonomy.new( "dog", "tax_id", "tax_source" )
                seqn = Sequence.new( "seqn", "VVVVV", "acc", "acc source", tax, "symbol", "2accession", "2source" )
                seqc = seqn.copy
                if ( seqc.get_name() != "seqn" )
                    return false
                end
                if ( seqc.get_accession() != "acc" )
                    return false
                end
                if ( seqc.get_accession_source() != "acc source" )
                    return false
                end
                if ( seqc.get_taxonomy.get_name != "dog" )
                    return false
                end
                if ( seqc.get_taxonomy.get_id != "tax_id" )
                    return false
                end
                if ( seqc.get_symbol != "symbol" )
                    return false
                end
                if ( seqc.get_secondary_accession != "2accession" )
                    return false
                end
                if ( seqc.get_secondary_accession_source != "2source" )
                    return false
                end

            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end


        def test_msa()
            begin
                msa = Msa.new()
                seq0 = Sequence.new( "seq 0", "a-*-_ x-ijklmnopqrstuvwxyz" )
                seq1 = Sequence.new( "seq 1", "ab--_ X-ijklmnopqrstuvwxyz" )
                seq2 = Sequence.new( "seq 2", "abc-_?x-ijklmnopqrstuvwxyz" )
                seq3 = Sequence.new( "seq 3", "abcd_?x-ijklmnopqrstuvwxyz" )
                seq4 = Sequence.new( "seq 4", "abcde?x-ijklmnopqrstuvwxyz" )
                seq5 = Sequence.new( "seq 5", "abcdefx-ijklmnopqrstuvwxyz" )
                msa.add_sequence( seq0 );
                msa.add_sequence( seq1 );
                msa.add_sequence( seq2 );
                msa.add_sequence( seq3 );
                msa.add_sequence( seq4 );
                msa.add_sequence( seq5 );
                msa.add(                      "seq 6", "abcdefg-ijklmnopqrstuvwxyz" );
                if ( msa.get_sequence( 0 ).get_name() != "seq 0" )
                    return false
                end
                if ( msa.get_by_name( "Eq 1", false, true ).get_name != "seq 1" )
                    return false
                end
                if ( msa.find_by_name( "Eq 2", false, true )[ 0 ] != 2 )
                    return false
                end
                if ( !msa.is_aligned )
                    return false
                end
                if ( msa.get_number_of_seqs != 7 )
                    return false
                end
                if ( msa.get_length != 26 )
                    return false
                end
                msa.add( "seq 7", "abcdefgqijklmnopqrstuvwxyz" );
                if ( msa.get_number_of_seqs != 8 )
                    return false
                end
                msa.remove_sequence!( 7 )
                if ( msa.get_number_of_seqs != 7 )
                    return false
                end
                msa.remove_gap_only_columns!()
                if ( msa.get_length() != 25 )
                    return false
                end
                if ( msa.get_by_name( "seq 0" ).get_sequence_as_string != "a-*-_ xijklmnopqrstuvwxyz" )
                    return false
                end
                msa.remove_gap_columns_w_gap_ratio!( 6.1 / 7.0 )
                if ( msa.get_length() != 25 )
                    return false
                end
                msa.remove_gap_columns_w_gap_ratio!( 6.0 / 7.0 )
                if ( msa.get_length() != 25 )
                    return false
                end
                if ( msa.get_by_name( "seq 0" ).get_sequence_as_string != "a-*-_ xijklmnopqrstuvwxyz" )
                    return false
                end
                msa.remove_gap_columns_w_gap_ratio!( 5.0 / 7.0 )
                if ( msa.get_length() != 25 )
                    return false
                end
                if ( msa.get_by_name( "seq 0" ).get_sequence_as_string != "a-*-_ xijklmnopqrstuvwxyz" )
                    return false
                end
                msa.remove_gap_columns_w_gap_ratio!( 2.0 / 7.0 )
                if ( msa.get_length() != 23 )
                    return false
                end
                if ( msa.get_by_name( "seq 0" ).get_sequence_as_string != "a-* xijklmnopqrstuvwxyz" )
                    puts msa.get_by_name( "seq 0" ).get_sequence_as_string
                    return false
                end
                msa.remove_gap_columns_w_gap_ratio!( 1.0 / 7.0 )
                if ( msa.get_length() != 21 )
                    return false
                end
                if ( msa.get_by_name( "seq 0" ).get_sequence_as_string != "a-xijklmnopqrstuvwxyz" )
                    return false
                end
                msa2 = Evoruby::Msa.new()
                msa2.add( "seq0", "abcdefgh" );
                msa2.add( "seq1", "a-cdefgh" );
                msa2.add( "seq2", "a--defgh" );
                msa2.add( "seq3", "a---efgh" );
                msa2.add( "seq4", "a----fgh" );
                msa2.add( "seq5", "a" );
                if ( msa2.is_aligned )
                    return false
                end
                msa2.remove_sequence!( 5 )
                if ( !msa2.is_aligned )
                    return false
                end
                if ( msa2.get_number_of_seqs != 5 )
                    return false
                end
                msa2.remove_gap_only_columns!()

                if ( msa2.get_length != 8 )
                    return false
                end

                msa2.remove_sequences_by_gap_ratio!( 4.0 / 8.0 )
                if ( msa2.get_number_of_seqs != 5 )
                    return false
                end
                msa2.remove_sequences_by_gap_ratio!( 3.0 / 8.0 )
                if ( msa2.get_number_of_seqs != 4 )
                    return false
                end
                msa2.remove_sequences_by_gap_ratio!( 1.0 / 8.0 )
                if ( msa2.get_number_of_seqs != 2 )
                    return false
                end
                msa2.remove_sequences_by_gap_ratio!( 0.0 )
                if ( msa2.get_number_of_seqs != 1 )
                    return false
                end
                msa2.add( "seq1", "a-cdefgh" );
                msa2.add( "seq2", "a--defgh" );
                msa2.add( "seq3", "a---efgh" );
                msa2.add( "seq4", "a----fgh" );

                msa2.remove_sequences_by_non_gap_length!( 4 )
                if ( msa2.get_number_of_seqs != 5 )
                    return false
                end
                msa2.remove_sequences_by_non_gap_length!( 5 )
                if ( msa2.get_number_of_seqs != 4 )
                    return false
                end
                msa2.remove_sequences_by_non_gap_length!( 8 )
                if ( msa2.get_number_of_seqs != 1 )
                    return false
                end
                msa2.add( "seq1", "a-cdefgh" );
                msa2.add( "seq2", "a--defgh" );
                msa2.add( "seq3", "a---efgh" );
                msa2.add( "seq4", "a----fgh" );
                msa2.trim!( 0, 7 )
                if ( msa2.get_by_name( "seq0" ).get_sequence_as_string != "abcdefgh" )
                    return false
                end
                msa2.trim!( 3, 4 )
                if ( msa2.get_by_name( "seq0" ).get_sequence_as_string != "de" )
                    return false
                end
                msa3 = Evoruby::Msa.new()
                msa3.add( "seq0", "abcdefgh-abcdef--*" );
                msa3.add( "seq1", "b-deefgh-a____f--*" );
                msa3.add( "seq2", "A________abcdef--*" );
                msa3.add( "seq3", "A   Efgh---------*" );
                msa3.add( "seq4", "    eFhh---------*" );
                msa3.add( "seq5", "----------------ee" );
                if ( !Test::same?( msa3.calculate_overlap( 0, 0 ), 14 ) )
                    return false
                end
                if ( !Test::same?( msa3.calculate_overlap( 0, 1 ), 9 ) )
                    return false
                end
                if ( !Test::same?( msa3.calculate_overlap( 0, 5 ), 0 ) )
                    return false
                end
                if ( !Test::same?( msa3.calculate_overlap( 4, 5 ), 0 ) )
                    return false
                end
                if ( !msa3.overlap?( 2, 3 ) )
                    return false
                end
                if ( msa3.overlap?( 2, 3, 2 ) )
                    return false
                end
                if ( msa3.overlap?( 4, 5 ) )
                    return false
                end
                if ( !Test::same?( msa3.calculate_identities( 4, 5 ), 0 ) )
                    return false
                end
                if ( !Test::same?( msa3.calculate_identities( 3, 4 ), 3 ) )
                    return false
                end
                if ( msa3.split_into_overlapping_msa.length != 3 )
                    return false
                end
                if ( msa3.split_into_overlapping_msa( 5 ).length != 4 )
                    return false
                end


                msa4 = Msa.new()
                seq0 = Sequence.new( "seq 0", "ABCDED" )
                seq1 = Sequence.new( "seq 1", "ABCDEE" )
                seq2 = Sequence.new( "seq 2", "abcded" )
                seq3 = Sequence.new( "seq 3", " ABCDEE" )
                seq4 = Sequence.new( "seq 4", "ABCDEV" )
                seq5 = Sequence.new( "seq 5", "ABCDED" )
                seq6 = Sequence.new( "seq 6", "AB.DEI" )
                seq7 = Sequence.new( "seq 7", "aB-DEi*" )
                seq8 = Sequence.new( "seq 8", "ABCDED" )
                seq9 = Sequence.new( "seq 9", "ABCDED" )
                seq10 = Sequence.new( "seq 10", "ABCDED" )
                seq11 = Sequence.new( "seq 11", "ABCDED" )
                msa4.add_sequence( seq0 );
                msa4.add_sequence( seq1 );
                msa4.add_sequence( seq2 );
                msa4.add_sequence( seq3 );
                msa4.add_sequence( seq4 );
                msa4.add_sequence( seq5 );
                msa4.add_sequence( seq6 );
                msa4.add_sequence( seq7 );
                msa4.add_sequence( seq8 );
                msa4.add_sequence( seq9 );
                msa4.add_sequence( seq10 );
                msa4.add_sequence( seq11 );

                msa4.remove_redundant_sequences!

                puts msa4.to_str

                if msa4.get_number_of_seqs != 4
                    return false
                end

                if msa4.get_sequence( 0 ).get_name != "seq 0"
                    return false
                end
                if msa4.get_sequence( 1 ).get_name != "seq 1"
                    return false
                end
                if msa4.get_sequence( 2 ).get_name != "seq 4"
                    return false
                end
                if msa4.get_sequence( 3 ).get_name != "seq 6"
                    return false
                end

            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_msa_factory()
            begin
                f = MsaFactory.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_domain_structure()
            begin
                ds = DomainStructure.new( 190 )
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_protein_domain()
            begin
                ds = ProteinDomain.new( "domain", 23, 466, "d1", 0.4 )
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_basic_table()
            begin
                t = BasicTable.new()
                t.set_value( 233, 923, "snake" )
                t.set_value( 233, 923, "lizard" )
                if ( t.get_value_as_string( 233, 923 ) != "lizard" )
                    return false
                end
                if ( t.get_value_as_string( 33, 23 ) != "" )
                    return false
                end
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_msa_io()
            begin
                msaio = MsaIO.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_phylip_sequentialwriter()
            begin
                p = PhylipSequentialWriter.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_nexus_writer()
            begin
                n = NexusWriter.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_fasta_writer()
            begin
                f = FastaWriter.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_general_msa_parser( path_to_evoruby )
            begin
                g = GeneralMsaParser.new()
                f = MsaFactory.new()
                sep = ""
                if ( !Util::is_string_empty?( path_to_evoruby ) )
                    sep = Constants::FILE_SEPARATOR
                end
                msa = f.create_msa_from_file( path_to_evoruby +
                     sep +
                     GENERAL_MSA_FILE, g )

                if ( msa.get_length() != 29 )
                    return false
                end
                if ( msa.get_number_of_seqs() != 7 )
                    return false
                end

                seq0 = msa.get_sequence( 0 )
                seq1 = msa.get_sequence( 1 )
                seq2 = msa.get_sequence( 2 )
                seq3 = msa.get_sequence( 3 )
                seq4 = msa.get_sequence( 4 )
                seq5 = msa.get_sequence( 5 )
                seq6 = msa.get_sequence( 6 )

                if ( seq0.get_name() != "sequence0" )
                    return false
                end
                if ( seq0.get_sequence_as_string() != "ABCDE.GHIJKLMNOPQR.TUVWabcxy0" )
                    return false
                end

                if ( seq1.get_name() != "sequence1" )
                    return false
                end
                if ( seq1.get_sequence_as_string() != "abcdefghijklmnopqrstuvwabcxy1" )
                    return false
                end

                if ( seq2.get_name() != "sequence2" )
                    return false
                end
                if ( seq2.get_sequence_as_string() != "abcdefghijkl---x_-*?_XXabcxy2" )
                    return false
                end

                if ( seq3.get_name() != "sequence3" )
                    return false
                end
                if ( seq3.get_sequence_as_string() != "12345678901234567890123abcxy3" )
                    return false
                end

                if ( seq4.get_name() != "sequence4" )
                    return false
                end
                if ( seq4.get_sequence_as_string() != "--------------------------xy4" )
                    return false
                end

                if ( seq5.get_name() != "sequence5" )
                    return false
                end
                if ( seq5.get_sequence_as_string() != "a*c*ef****************wabcxy5" )
                    return false
                end

                if ( seq6.get_name() != "sequence6" )
                    return false
                end
                if ( seq6.get_sequence_as_string() != "ururufhfghfgftgfhftgfttabcxy6" )
                    return false
                end

            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_fasta_parser( path_to_evoruby )
            begin
                fasta = FastaParser.new()
                f = MsaFactory.new()
                sep = ""
                if ( !Util::is_string_empty?( path_to_evoruby ) )
                    sep = Constants::FILE_SEPARATOR
                end
                msa = f.create_msa_from_file( path_to_evoruby +
                     sep +
                     FASTA_FILE, fasta )

                if ( msa.get_length() != 6 )
                    return false
                end
                if ( msa.get_number_of_seqs() != 4 )
                    return false
                end

                seq0 = msa.get_sequence( 0 )
                seq1 = msa.get_sequence( 1 )
                seq2 = msa.get_sequence( 2 )
                seq3 = msa.get_sequence( 3 )

                if ( seq0.get_name() != "sequence 0" )
                    return false
                end
                if ( seq0.get_sequence_as_string() != "ABCDEF" )
                    return false
                end

                if ( seq1.get_name() != "sequence 1" )
                    return false
                end
                if ( seq1.get_sequence_as_string() != "abcdef" )
                    return false
                end

                if ( seq2.get_name() != "sequence 2" )
                    return false
                end
                if ( seq2.get_sequence_as_string() != "123456" )
                    return false
                end
                if ( seq3.get_name() != "sequence 3" )
                    return false
                end
                if ( seq3.get_sequence_as_string() != "a-c--f" )
                    return false
                end

            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_ncbi_tseq_parser( path_to_evoruby )
            begin
                parser = NcbiTSeqParser.new
                f = MsaFactory.new
                sep = ""
                if ( !Util::is_string_empty?( path_to_evoruby ) )
                    sep = Constants::FILE_SEPARATOR
                end
                msa = f.create_msa_from_file( path_to_evoruby +
                     sep +
                     TSEQ_FILE, parser )

                if ( msa.get_number_of_seqs() != 9 )
                    return false
                end

                seq0 = msa.get_sequence( 0 )
                seq1 = msa.get_sequence( 1 )
                seq8 = msa.get_sequence( 8 )

                if ( seq0.get_name() != "SusD [Bacteroides thetaiotaomicron VPI-5482]" )
                    return false
                end
                if ( seq0.get_sequence_as_string() != "MKTKYIKQLFSAALIAVLSSGVTSCINDLDISPIDPQTGGSFDQQGVFVKGYAMLGVTGQKGIDGSPDLDGQDEGESGFYRTTFNCNELPTDECLWAWQENQDIPQLTSISWSPSSQRTEWVYVRLGYDITQYNFFLDQTEGMTDAETLRQRAEIRFLRALHYWYFLDLFGKAPFKEHFSNDLPVEKKGTELYTYIQNELNEIEADMYEPRQAPFGRADKAANWLLRARLYLNAGVYTGQTDYAKAEEYASKVIGSAYKLCTNYSELFMADNDENENAMQEIILPIRQDGVKTRNYGGSTYLVCGTRVAGMPRMGTTNGWSCIFARAAMVQKFFSNLEDVPMLPADVEIPTKGLDTDEQIDAFDAEHGIRTEDMIKAAGDDRALLYSGVGGGRRKIQTDAISGFTDGLSIVKWQNYRSDGKPVSHATYPDTDIPLFRLAEAYLTRAEAIFRQGGDATGDINELRKRANCTRKVQTVTEQELIDEWAREFYLEGRRRSDLVRFGMFTTNKYLWDWKGGAMNGTSVASYYNKYPIPVSDINNNRNMSQNEGYK" )
                    return false
                end
                if ( seq0.get_accession != "29341016" )
                    return false
                end
                if ( seq0.get_accession_source != "gi" )
                    return false
                end
                if ( seq0.get_taxonomy.get_name != "Bacteroides thetaiotaomicron VPI-5482" )
                    return false
                end
                if ( seq0.get_taxonomy.get_id != "226186" )
                    return false
                end
                if ( seq0.get_taxonomy.get_id_source != "ncbi" )
                    return false
                end


                if ( seq1.get_name() != "SusD, outer membrane protein [Bacteroides thetaiotaomicron VPI-5482]" )
                    return false
                end
                if ( seq1.get_accession != "29349109" )
                    return false
                end
                if ( seq1.get_accession_source != "gi" )
                    return false
                end
                if ( seq1.get_taxonomy.get_name != "Bacteroides thetaiotaomicron VPI-5482" )
                    return false
                end
                if ( seq1.get_taxonomy.get_id != "226186" )
                    return false
                end
                if ( seq1.get_taxonomy.get_id_source != "ncbi" )
                    return false
                end


                if ( seq8.get_name() != "Chain A, B. Thetaiotaomicron Susd With Maltotriose" )
                    return false
                end
                if ( seq8.get_accession != "pdb|3CKB|A" )
                    return false
                end
                if ( seq8.get_accession_source != "ncbi" )
                    return false
                end
                if ( seq8.get_taxonomy.get_name != "Bacteroides thetaiotaomicron" )
                    return false
                end
                if ( seq8.get_taxonomy.get_id != "818" )
                    return false
                end
                if ( seq8.get_taxonomy.get_id_source != "ncbi" )
                    return false
                end

            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_hmmsearch_domain_extractor()
            begin
                h = Evoruby::HmmsearchDomainExtractor.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_domain_sequence_extractor()
            begin
                h = Evoruby::DomainSequenceExtractor.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_hmmscan_parser()
            begin
                h = Evoruby::HmmscanSummary.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_domains_to_forester()
            begin
                d = Evoruby::DomainsToForester.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end


        def test_basic_table_parser()
            begin
                b = Evoruby::BasicTableParser.new()
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end


        def test_cla()
            begin
                cla = CommandLineArguments.new( Array.new )
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_tree_puzzle()
            begin
                tp = TreePuzzle.new()
                tp.run( '/home/czmasek/scratch/small.aln',
                    :wag,
                    :uniform,
                    200 )
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end

        def test_fastme()
            begin
                fastme = FastMe.new()
                fastme.run( '/home/czmasek/scratch/outdist', 0, :GME )
            rescue Exception => e
                puts()
                puts( e.to_s )
                puts()
                return false
            end
            return true
        end


        def run()

            t0 = Time.now

            puts
            puts "ruby version " + RUBY_VERSION
            puts Constants::EVORUBY + " version " + Constants::EVORUBY_VERSION
            puts

            path_to_evoruby = Test.get_path_to_evoruby()

            if ( Util.is_string_empty?( path_to_evoruby ) )
                path_to_evoruby = ""
                puts()
                puts( "Warning! Path to evoruby could not be established. Some tests will might fail." )
                puts()
            end

            print( "--- Taxonomy: " )
            if ( test_taxonomy() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- Sequence: " )
            if ( test_sequence() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- Msa: " )
            if ( test_msa() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- MsaFactory: " )
            if ( test_msa_factory() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- DomainStructure: " )
            if ( test_domain_structure() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- ProteinDomain: " )
            if ( test_protein_domain() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- BasicTable: " )
            if ( test_basic_table() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- MsaIO: " )
            if ( test_msa_io )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- PhylipSequentialWriter: " )
            if ( test_phylip_sequentialwriter )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- FastaWriter : " )
            if ( test_fasta_writer )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- NexusWriter: " )
            if ( test_nexus_writer )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- FastaParser: " )
            if ( test_fasta_parser( path_to_evoruby ) )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- NCBI Tseq parser: " )
            if (  test_ncbi_tseq_parser( path_to_evoruby ) )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- GeneralMsaParser: " )
            if ( test_general_msa_parser( path_to_evoruby ) )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end


            print( "--- Hmmsearch domain extractor: " )
            if ( test_hmmsearch_domain_extractor )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- Domain sequence extractor: " )
            if ( test_domain_sequence_extractor )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- Hmmscan parser: " )
            if ( test_hmmscan_parser )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end


            print( "--- Domains 2 forester: " )
            if ( test_domains_to_forester )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- BasicTableParser: " )
            if ( test_basic_table_parser )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- TreePuzzle (wrapper): " )
            if ( test_tree_puzzle() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end

            print( "--- FastMe (wrapper): " )
            if ( test_fastme() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end



            print( "--- CLA: " )
            if ( test_cla() )
                puts( "ok" )
                @successes += 1
            else
                puts( "FAILED" )
                @failures += 1
            end
            puts
            puts "ruby version " + RUBY_VERSION
            puts Constants::EVORUBY + " version " + Constants::EVORUBY_VERSION
            puts

            td = Time.at( Time.now - t0 )
            puts( "Time            : #{ td.sec }.#{ td.usec }s" )
            puts()

            puts( "Successful tests: " + @successes.to_s )
            puts( "Failed tests    : " + @failures.to_s )
            puts()
            if ( @failures < 1 )
                puts( "OK" )
            else
                puts( "NOT ok" )
            end

            puts()
        end

        private

        def Test.same?( n, m )
            return ( ( n - m ).abs < 0.000001 )
        end

        def Test.get_path_to_evoruby()
            rubylib = ENV['RUBYLIB'].split(':')
            evoruby_path = nil
            rubylib.each do | path |
                if ( path =~ /evoruby/ )
                    evoruby_path = path
                    break
                end
            end
            evoruby_path
        end

    end # class Test


    test = Test.new()

    test.run()


end # module Evoruby


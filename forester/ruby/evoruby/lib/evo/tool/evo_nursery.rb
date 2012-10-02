#
# = lib/evo/apps/evo_nursery.rb - EvoNursery class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: evo_nursery.rb,v 1.11 2010/12/13 19:00:11 cmzmasek Exp $



require 'lib/evo/soft/fastme'
require 'lib/evo/soft/tree_puzzle'
require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/writer/phylip_sequential_writer'
require 'lib/evo/io/parser/general_msa_parser'
require 'lib/evo/io/writer/msa_writer'

require 'iconv'

module Evoruby

    class EvoNursery
        GAP_RATIO           = 0.50
        GAP_RATIO_FOR_SEQS  = 0.75
        MIN_LENGTH          = 30
        MIN_SEQS            = 4
        MAX_SEQS            = 3000
        MAX_ALN_FILE_SIZE   = 5000000
        MODEL               = :auto
        RATES               = :uniform
        FASTME_INITIAL_TREE = :GME
        ALN_NAME            = '_align_'
        TREE_PUZZLE_OUTDIST = TreePuzzle::OUTDIST
        TREE_PUZZLE_OUTFILE = TreePuzzle::OUTFILE
        FASTME_OUTTREE      = FastMe::OUTTREE
        FASTME_OUTPUT_D     = FastMe::OUTPUT_D

        PRG_NAME       = "evo_nursery"
        PRG_DATE       = "2012.03.21"
        PRG_DESC       = "pfam alignments to evolutionary trees"
        PRG_VERSION    = "0.20"
        COPYRIGHT      = "2009-2012 Christian M Zmasek"
        CONTACT        = "phylosoft@gmail.com"
        WWW            = "www.phylosoft.org"

        HELP_OPTION_1       = "help"
        HELP_OPTION_2       = "h"

        def run

            Util.print_program_information( PRG_NAME,
                PRG_VERSION,
                PRG_DESC,
                PRG_DATE,
                COPYRIGHT,
                CONTACT,
                WWW,
                STDOUT )

            if RUBY_VERSION !~ /1.9/
                puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
                exit( -1 )
            end

            forester_home = Util.get_env_variable_value( Constants::FORESTER_HOME_ENV_VARIABLE )
            java_home = Util.get_env_variable_value( Constants::JAVA_HOME_ENV_VARIABLE )
            decorator = java_home + '/bin/java -cp ' + forester_home + '/java/forester.jar org.forester.application.decorator'

            if ( ARGV == nil || ARGV.length != 1 )
                help
                exit( -1 )
            end

            begin
                cla = CommandLineArguments.new( ARGV )
            rescue ArgumentError => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s )
            end

            if ( cla.is_option_set?( HELP_OPTION_1 ) ||
                     cla.is_option_set?( HELP_OPTION_2 ) )
                help
                exit( 0 )
            end

            output_dir = cla.get_file_name( 0 )

            if output_dir !~ /\/$/
                output_dir = output_dir + '/'
            end

            if !File.exists?( output_dir )
                Util.fatal_error( PRG_NAME, output_dir.to_s + " does not exist", STDOUT )
            end
            ic = Iconv.new( 'UTF-8//IGNORE', 'UTF-8' )
            files = Dir.entries( "." )
            skipped = Array.new
            counter = 1
            analyzed = 0;
            begin
                files.each { |pfam_aln_file|
                    if ( !File.directory?( pfam_aln_file ) &&
                             pfam_aln_file !~ /^\./ &&
                             pfam_aln_file !~ /.+\.tre$/  )

                        tree_out_file = output_dir + File.basename( pfam_aln_file ) + ".xml"

                        if File.exists?( tree_out_file )
                            puts
                            puts
                            puts "***** skipping " + File.basename( pfam_aln_file ) + ", already exists"
                            puts
                            skipped.push( File.basename( pfam_aln_file ) + " [already exists]" )
                            next
                        end

                        puts
                        puts counter.to_s + ": " + pfam_aln_file.to_str
                        counter += 1
                        if File.size( pfam_aln_file ) > MAX_ALN_FILE_SIZE
                            puts "***** skipping, file size: " +  File.size( pfam_aln_file ).to_s
                            skipped.push( File.basename( pfam_aln_file ) + " [file size: " +  File.size( pfam_aln_file ).to_s + "]" )
                            next
                        end

                        f = MsaFactory.new()
                        msa = f.create_msa_from_file( pfam_aln_file, GeneralMsaParser.new() )

                        if msa.get_number_of_seqs < MIN_SEQS || msa.get_number_of_seqs > MAX_SEQS
                            puts "***** skipping, seqs: " + msa.get_number_of_seqs.to_s
                            skipped.push( File.basename( pfam_aln_file ) + " [seqs: " +  msa.get_number_of_seqs.to_s + "]" )
                            next
                        end

                        msa.remove_gap_columns_w_gap_ratio!( GAP_RATIO )

                        length = msa.get_length
                        if length < MIN_LENGTH
                            puts "***** skipping, length: " + length.to_s
                            skipped.push( File.basename( pfam_aln_file ) + " [length: " +  length.to_s + "]" )
                            next
                        end

                        msa.remove_sequences_by_gap_ratio!( GAP_RATIO_FOR_SEQS )

                        if msa.get_number_of_seqs < MIN_SEQS
                            puts "***** skipping, seqs: " + msa.get_number_of_seqs.to_s
                            skipped.push( File.basename( pfam_aln_file ) + " [seqs: " +  msa.get_number_of_seqs.to_s + "]" )
                            next
                        end

                        map_file = output_dir + File.basename( pfam_aln_file ) + ".map"
                        f = File.open( map_file, 'a' )
                        for i in 0 ... msa.get_number_of_seqs
                            name = msa.get_sequence( i ).get_name()
                            name =~ /(.+)_(.+)\/.+/
                            acc = $1
                            tax_code = $2

                            mapping_str = i.to_s
                            mapping_str << "\t"
                            mapping_str << 'TAXONOMY_CODE:'
                            mapping_str << tax_code
                            mapping_str << "\t"
                            mapping_str << 'SEQ_SYMBOL:'
                            mapping_str << ( acc + '_' + tax_code )
                            mapping_str << "\t"
                            if ( acc.length < 6 )
                                acc = acc + '_' + tax_code
                            end
                            mapping_str << 'SEQ_ACCESSION:'
                            mapping_str << acc
                            mapping_str << "\t"
                            mapping_str << 'SEQ_ACCESSION_SOURCE:UniProtKB'
                            mapping_str << "\t"
                            mapping_str << 'NODE_NAME:'
                            mapping_str << name
                            f.print( mapping_str )
                            f.print( "\n" )
                            name = msa.get_sequence( i ).set_name( i.to_s )
                        end
                        f.close

                        io = MsaIO.new()
                        w = MsaWriter
                        w = PhylipSequentialWriter.new()
                        w.clean( true )
                        w.set_max_name_length( 10 )
                        if File.exists?( output_dir + ALN_NAME )
                            File.unlink( output_dir + ALN_NAME )
                        end
                        io.write_to_file( msa, output_dir + ALN_NAME, w )

                        tp = TreePuzzle.new()
                        tp.run( output_dir + ALN_NAME,
                            MODEL,
                            RATES,
                            msa.get_number_of_seqs )

                        File.rename( output_dir + ALN_NAME, output_dir  + File.basename( pfam_aln_file ) + ".aln" )

                        fastme = FastMe.new()
                        fastme.run( TREE_PUZZLE_OUTDIST, 0, FASTME_INITIAL_TREE )

                        pfam_acc = nil
                        pfam_de = nil
                        File.open( pfam_aln_file ) do |file|
                            while line = file.gets
                                line = ic.iconv( line )
                                if line =~ /^#=AC\s+(.+)/
                                    pfam_acc = $1
                                end
                                if line =~ /^#=DE\s+(.+)/
                                    pfam_de = $1
                                end
                                if pfam_acc && pfam_de
                                    break
                                end
                            end
                        end
                        if !pfam_acc || !pfam_de
                            Util.fatal_error( PRG_NAME, "problem with " + pfam_aln_file.to_s, STDOUT )
                        end

                        puzzle_model = nil
                        File.open( TREE_PUZZLE_OUTFILE ) do |file|
                            while line = file.gets
                                line = ic.iconv( line )
                                if line =~ /^Model\s+of\s+substitution:\s+(.+)/
                                    puzzle_model = $1
                                    break
                                end
                            end
                        end
                        if !puzzle_model
                            Util.fatal_error( PRG_NAME, "problem with puzzle outfile: " + TREE_PUZZLE_OUTFILE.to_s, STDOUT )
                        end

                        desc = pfam_de
                        desc << ' | '
                        desc << 'ML pwd estimation by TREE-PUZZLE version '
                        desc << TreePuzzle::VERSION
                        desc << ', model: '
                        desc << puzzle_model
                        desc << ', rates: '
                        desc << RATES.to_s
                        desc << '; tree estimation by FastME version '
                        desc << FastMe::VERSION
                        desc << ', initial tree: '
                        desc << FASTME_INITIAL_TREE.to_s
                        desc << '; aln length: '
                        desc << msa.get_length.to_s

                        cmd = decorator + " -table -p -pn=\"" + pfam_aln_file +
                         "\" -pi=pfam:" + pfam_acc +
                         " -pd=\"" + desc + "\" " +
                         FASTME_OUTTREE + ' ' +
                         map_file + ' ' + tree_out_file

                        IO.popen( cmd , 'r+' ) do | pipe |
                            pipe.close_write
                        end
                        analyzed += 1

                        File.unlink( map_file )
                        File.unlink(TREE_PUZZLE_OUTDIST)
                        File.unlink( TREE_PUZZLE_OUTFILE )
                        File.unlink( FASTME_OUTPUT_D )
                    end
                }
            rescue ArgumentError, IOError, StandardError => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
            end

            puts()
            puts( 'Skipped:' )
            puts()
            for i in 0 ... skipped.size
                puts i.to_s + ": " + skipped[ i ]
            end

            puts()
            puts( 'Skipped : ' + skipped.size.to_s + ' alignments' )
            puts( 'Analyzed: ' +  analyzed.to_s    + ' alignments' )

            puts( 'Min gap ratio for col del  : ' + GAP_RATIO.to_s )
            puts( 'Min gap ratio for seq del  : ' + GAP_RATIO_FOR_SEQS.to_s )
            puts( 'Minimal aln length         : ' + MIN_LENGTH.to_s )
            puts( 'Minimal number of sequences: ' + MIN_SEQS.to_s )
            puts( 'Maximal number of sequences: ' + MAX_SEQS.to_s )
            puts( 'Maximal aln file size      : ' + MAX_ALN_FILE_SIZE.to_s )
            puts( 'Model              : ' + MODEL.to_s )
            puts( 'FastME initial tree: ' + FASTME_INITIAL_TREE.to_s )

            puts()
            puts( '[' + PRG_NAME + '] > OK' )
            puts()

        end  # run

        private

        def help
            puts( "Usage:" )
            puts()
            puts( "  " + PRG_NAME + ".rb <output dir> " )
            puts()
        end


    end # class EvoNursery

end # module Evoruby
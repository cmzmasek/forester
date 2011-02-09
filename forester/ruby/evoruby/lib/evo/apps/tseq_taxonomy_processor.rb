#
# = lib/evo/apps/tseq_taxonomy_processor - TseqTaxonomyProcessor class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: tseq_taxonomy_processor.rb,v 1.6 2010/12/13 19:00:11 cmzmasek Exp $


require 'lib/evo/util/util'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/msa/msa'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/parser/sp_taxonomy_parser'
require 'lib/evo/io/parser/ncbi_tseq_parser'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/writer/phylip_sequential_writer'
require 'lib/evo/util/command_line_arguments'

module Evoruby

    class TseqTaxonomyProcessor

        PRG_NAME       = "tseq_tap"
        PRG_DATE       = "2009.01.06"
        PRG_DESC       = "preprocessing of multiple sequence files in ncbi tseq xml format"
        PRG_VERSION    = "1.02"
        COPYRIGHT      = "2009 Christian M Zmasek"
        CONTACT        = "phylosoft@gmail.com"
        WWW            = "www.phylosoft.org"

        TAXONOMY_CODE           = "TAXONOMY_CODE:"
        TAXONOMY_ID             = "TAXONOMY_ID:"
        TAXONOMY_ID_TYPE        = "TAXONOMY_ID_TYPE:"
        TAXONOMY_SN             = "TAXONOMY_SN:"
        TAXONOMY_CN             = "TAXONOMY_CN:"
        SEQ_ACCESSION           = "SEQ_ACCESSION:"
        SEQ_ACCESSION_SOURCE    = "SEQ_ACCESSION_SOURCE:"
        SEQ_SECONDARY_ACCESSION = "SEQ_SECONDARY_ACCESSION:"
        SEQ_SYMBOL              = "SEQ_SYMBOL:"
        SEQ_NAME                = "SEQ_NAME:"
        SEQ_MOL_SEQ             = "SEQ_MOL_SEQ:"

        def initialize()
            @tax_ids_to_sp_taxonomies = Hash.new()
        end

        def run()

            Util.print_program_information( PRG_NAME,
                PRG_VERSION,
                PRG_DESC,
                PRG_DATE,
                COPYRIGHT,
                CONTACT,
                WWW,
                STDOUT )

            if  ARGV == nil || ARGV.length != 4
                puts( "Usage: #{PRG_NAME}.rb <sp taxonomy file> <sequences in tseq xml format> <name for fasta outfile> <name for map outfile>" )
                puts()

                exit( -1 )
            end

            begin
                cla = CommandLineArguments.new( ARGV )
            rescue ArgumentError => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s )
            end
            allowed_opts = Array.new
            disallowed = cla.validate_allowed_options_as_str( allowed_opts )
            if ( disallowed.length > 0 )
                Util.fatal_error( PRG_NAME, "unknown option(s): " + disallowed )
            end

            sp_taxonomy_infile = cla.get_file_name( 0 )
            sequences_infile = cla.get_file_name( 1 )
            sequences_outfile = cla.get_file_name( 2 )
            mapping_outfile = cla.get_file_name( 3 )

            Util.fatal_error_if_not_readable( PRG_NAME, sp_taxonomy_infile )
            Util.fatal_error_if_not_readable( PRG_NAME, sequences_infile )
            Util.fatal_error_if_not_writable( PRG_NAME, mapping_outfile )
            Util.fatal_error_if_not_writable( PRG_NAME, sequences_outfile )

            sp_taxonomies = SpTaxonomyParser.parse( sp_taxonomy_infile )

            Util.print_message( PRG_NAME, "read in taxonomic data for " + sp_taxonomies.size.to_s + " species from: " + sp_taxonomy_infile )

            tseq_parser = NcbiTSeqParser.new
            msa_fac = MsaFactory.new

            seqs = msa_fac.create_msa_from_file( sequences_infile, tseq_parser )

            Util.print_message( PRG_NAME, "read in " + seqs.get_number_of_seqs.to_s + " sequences from: " + sequences_infile )

            removed = seqs.remove_redundant_sequences!( true, true )

            if removed.size > 0
                Util.print_message( PRG_NAME, "going to ignore the following " + removed.size.to_s + " redundant sequences:" )
                removed.each { | seq_name |
                    puts seq_name
                }
                Util.print_message( PRG_NAME, "will process " + seqs.get_number_of_seqs.to_s + " non-redundant sequences" )
            end

            mapping_out = File.open( mapping_outfile, "a" )

            for i in 0 ... seqs.get_number_of_seqs
                seq = seqs.get_sequence( i )
                if seq.get_taxonomy == nil
                    Util.fatal_error( PRG_NAME, "sequence [" + seq.get_name + "] has no taxonomy information" )
                end
                seq.set_name( Util::normalize_seq_name( modify_name( seq, i, sp_taxonomies, mapping_out ), 10 ) )
            end

            io = MsaIO.new()

            w = FastaWriter.new()

            w.set_max_name_length( 10 )
            w.clean( true )
            begin
                io.write_to_file( seqs, sequences_outfile, w )
            rescue Exception => e
                Util.fatal_error( PRG_NAME, "failed to write file: " + e.to_s )
            end
            mapping_out.close()

            Util.print_message( PRG_NAME, "wrote: " + mapping_outfile )
            Util.print_message( PRG_NAME, "wrote: " + sequences_outfile )
            Util.print_message( PRG_NAME, "OK" )

        end

        private

        def modify_name( seq, i, sp_taxonomies, mapping_outfile )

            tax_id = seq.get_taxonomy.get_id
            matching_sp_taxonomy = nil

            if @tax_ids_to_sp_taxonomies.has_key?( tax_id )
                # This is so that a second lookup will be much faster.
                matching_sp_taxonomy = @tax_ids_to_sp_taxonomies[ tax_id ]
            else
                sp_taxonomies.each { |sp_taxonomy|
                    if ( sp_taxonomy.id == tax_id )
                        if  matching_sp_taxonomy != nil
                            Util.fatal_error( PRG_NAME, "taxonomy id [" + tax_id.to_s + "] is not unique" )
                        end
                        matching_sp_taxonomy = sp_taxonomy
                        @tax_ids_to_sp_taxonomies[ tax_id ] = sp_taxonomy
                    end
                }
            end
            if  matching_sp_taxonomy == nil
                Util.fatal_error( PRG_NAME, "taxonomy id [" + tax_id.to_s + "] for [" +  seq.get_taxonomy.get_name + "] not found" )
            end

            new_name = i.to_s( 16 ) + "_" + matching_sp_taxonomy.code

            seq_name = seq.get_name
            if  seq_name =~ /\[.+\]$/
                # Redundant taxonomy information hides here.
                seq_name = seq_name.sub(/\[.+\]$/, '')
            end
            if  seq_name =~ /^\s*hypothetical\s+protein\s*/i
                # Pointless information.
                seq_name = seq_name.sub( /^\s*hypothetical\s+protein\s*/i, '' )
            end

            mapping_outfile.print( new_name + "\t" +
                 TAXONOMY_CODE + matching_sp_taxonomy.code + "\t" +
                 TAXONOMY_ID + tax_id + "\t" +
                 TAXONOMY_ID_TYPE + seq.get_taxonomy.get_id_source + "\t" +
                 TAXONOMY_SN + matching_sp_taxonomy.scientific_name + "\t" +
                 SEQ_ACCESSION + seq.get_accession + "\t" +
                 SEQ_ACCESSION_SOURCE + seq.get_accession_source + "\t" +
                 SEQ_SYMBOL + seq.get_symbol + "\t" +
                 SEQ_NAME + seq_name + "\t" +
                 SEQ_MOL_SEQ + seq.get_sequence_as_string +
                 Constants::LINE_DELIMITER )
            new_name
        end

    end 

end # module Evoruby
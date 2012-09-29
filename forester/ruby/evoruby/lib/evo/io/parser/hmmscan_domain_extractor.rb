#
# = lib/evo/io/parser/hmmscan_domain_extractor.rb - HmmscanDomainExtractor class
#
# Copyright::  Copyright (C) 2012 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id:  $


require 'lib/evo/util/constants'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/parser/fasta_parser'


module Evoruby

  class HmmscanDomainExtractor

    TRIM_BY = 2

    def initialize
    end

    # raises ArgumentError, IOError, StandardError
    def parse( domain_id,
        hmmsearch_output,
        fasta_sequence_file,
        outfile,
        passed_seqs_outfile,
        failed_seqs_outfile,
        e_value_threshold,
        length_threshold,
        add_position,
        add_domain_number,
        add_domain_number_as_digit,
        add_domain_number_as_letter,
        trim_name,
        add_species,
        min_linker,
        log )

      Util.check_file_for_readability( hmmsearch_output )
      Util.check_file_for_readability( fasta_sequence_file )
      Util.check_file_for_writability( outfile )
      Util.check_file_for_writability( passed_seqs_outfile )
      Util.check_file_for_writability( failed_seqs_outfile )

      in_msa = nil
      factory = MsaFactory.new()
      in_msa = factory.create_msa_from_file( fasta_sequence_file, FastaParser.new() )

      if ( in_msa == nil || in_msa.get_number_of_seqs() < 1 )
        error_msg = "could not find fasta sequences in " + fasta_sequence_file
        raise IOError, error_msg
      end

      out_msa = Msa.new

      failed_seqs = Msa.new
      passed_seqs = Msa.new
      out_msa_pairs = nil
      out_msa_distance_partners = nil
      out_msa_singlets = nil
      if min_linker
        out_msa_pairs = Msa.new
        out_msa_distant_partners = Msa.new
        out_msa_singlets = Msa.new
      end

      ld = Constants::LINE_DELIMITER

      domain_pass_counter     = 0
      domain_fail_counter     = 0
      proteins_with_passing_domains = 0
      proteins_with_failing_domains = 0
      max_domain_copy_number_per_protein = -1
      max_domain_copy_number_sequence    = ""

      prev_sequence = nil
      prev_number   = nil
      prev_env_from = nil
      prev_env_to   = nil
      prev_i_e_value  = nil
      prev_is_pair = false

      File.open( hmmsearch_output ) do | file |
        while line = file.gets
          if !is_ignorable?( line ) && line =~ /^\S+\s+/

            #         tn      acc     tlen    query   acc     qlen    Evalue  score   bias    #       of      c-E     i-E     score   bias    hf      ht      af      at      ef      et      acc     desc
            #         1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21      22      23
            line =~ /^(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(.*)/

            target_name = $1
            if domain_id != target_name
              next
            end

            sequence = $4
            number   = $10.to_i
            out_of   = $11.to_i
            env_from = $20.to_i
            env_to   = $21.to_i
            i_e_value  = $13.to_f
            if ( number > max_domain_copy_number_per_protein )
              max_domain_copy_number_sequence    = sequence
              max_domain_copy_number_per_protein = number
            end
            if ( ( ( e_value_threshold < 0.0 ) || ( i_e_value <= e_value_threshold ) ) &&
                 ( ( length_threshold <= 0 )   || ( env_to - env_from + 1 ) >= length_threshold.to_f )  )

              extract_domain( sequence,
                number,
                out_of,
                env_from,
                env_to,
                in_msa,
                out_msa,
                add_position,
                add_domain_number,
                add_domain_number_as_digit,
                add_domain_number_as_letter,
                trim_name ,
                add_species )
              domain_pass_counter += 1

              if passed_seqs.find_by_name_start( sequence, true ).length < 1
                add_sequence( sequence, in_msa, passed_seqs )
                proteins_with_passing_domains += 1
              end

              if min_linker
                if ( ( e_value_threshold < 0.0 ) || ( prev_i_e_value <= e_value_threshold  ) ) &&
                   ( ( length_threshold <= 0 )   || (  ( prev_env_to - prev_env_from + 1 ) >= length_threshold.to_f    ) )

                  if sequence != prev_sequence
                    prev_is_pair = false
                  end

                  if out_of == 1

                    if sequence == prev_sequence
                      puts "sequence == prev_sequence && out_of == 1"
                      exit
                    end
                    extract_domain( sequence,
                      number,
                      out_of,
                      env_from,
                      env_to,
                      in_msa,
                      out_msa_singlets,
                      false,
                      true,
                      false,
                      false,
                      trim_name ,
                      add_species )

                  elsif sequence == prev_sequence

                    if  ( env_from - prev_env_to ) <= min_linker  #######
                      extract_domain( sequence,
                        prev_number.to_s + "+" + number.to_s,
                        out_of,
                        prev_env_from,
                        env_to,
                        in_msa,
                        out_msa_pairs,
                        false,
                        true,
                        false,
                        false,
                        trim_name ,
                        add_species )
                      prev_is_pair = true
                    else               #######
                      if !prev_is_pair
                        extract_domain( sequence,
                          prev_number,
                          out_of,
                          prev_env_from,
                          prev_env_to,
                          in_msa,
                          out_msa_distant_partners,
                          false,
                          true,
                          false,
                          false,
                          trim_name ,
                          add_species )
                      end
                      if number == out_of
                        extract_domain( sequence,
                          number,
                          out_of,
                          env_from,
                          env_to,
                          in_msa,
                          out_msa_distant_partners,
                          false,
                          true,
                          false,
                          false,
                          trim_name ,
                          add_species )
                      end
                      prev_is_pair = false
                    end                #######

                  end
                  prev_sequence = sequence
                  prev_number   = number
                  prev_env_from = env_from
                  prev_env_to   = env_to
                  prev_i_e_value  = i_e_value
                end

              else
                print( domain_fail_counter.to_s + ": " + sequence.to_s + " did not meet threshold(s)" )
                log << domain_fail_counter.to_s + ": " + sequence.to_s + " did not meet threshold(s)"
                if ( ( e_value_threshold.to_f >= 0.0 ) && ( i_e_value > e_value_threshold ) )
                  print( " iE=" + i_e_value.to_s )
                  log << " iE=" + i_e_value.to_s
                end
                if ( ( length_threshold.to_f > 0 ) && ( env_to - env_from + 1 ) < length_threshold.to_f )
                  le = env_to - env_from + 1
                  print( " l=" + le.to_s )
                  log << " l=" + le.to_s
                end
                print( Constants::LINE_DELIMITER )
                log << Constants::LINE_DELIMITER
                domain_fail_counter  += 1

                if failed_seqs.find_by_name_start( sequence, true ).length < 1
                  add_sequence( sequence, in_msa, failed_seqs )
                  proteins_with_failing_domains += 1
                end
              end
            end
          end # if !is_ignorable?( line ) && line =~ /^\S+\s+/
        end #  while line = file.gets
      end #   File.open( hmmsearch_output ) do | file |

      if domain_pass_counter < 1
        error_msg = "no domain sequences were extracted"
        raise StandardError, error_msg
      end

      log << Constants::LINE_DELIMITER
      puts( "Max domain copy number per protein : " + max_domain_copy_number_per_protein.to_s )
      log << "Max domain copy number per protein : " + max_domain_copy_number_per_protein.to_s
      log << Constants::LINE_DELIMITER

      if ( max_domain_copy_number_per_protein > 1 )
        puts( "First protein with this copy number: " + max_domain_copy_number_sequence )
        log << "First protein with this copy number: " + max_domain_copy_number_sequence
        log << Constants::LINE_DELIMITER
      end

      write_msa( out_msa, outfile  )
      write_msa( passed_seqs, passed_seqs_outfile )
      write_msa( failed_seqs, failed_seqs_outfile )

      if out_msa_pairs
        write_msa( out_msa_pairs, outfile +"_" + min_linker.to_s )
      end

      if out_msa_singlets
        write_msa( out_msa_singlets, outfile +"_singles" )
      end

      if out_msa_distant_partners
        write_msa( out_msa_distant_partners, outfile +"_singles" )
      end


      log << ld
      log << "passing domains              : " + domain_pass_counter.to_s + ld
      log << "failing domains              : " + domain_fail_counter.to_s + ld
      log << "proteins with passing domains: " + proteins_with_passing_domains.to_s + ld
      log << "proteins with failing domains: " + proteins_with_failing_domains.to_s + ld
      log << ld

      return domain_pass_counter

    end # parse


    private

    def write_msa( msa, filename )
      io = MsaIO.new()
      w = FastaWriter.new()
      w.set_line_width( 60 )
      w.clean( true )
      begin
        io.write_to_file( msa, filename, w )
      rescue Exception
        error_msg = "could not write to \"" + filename + "\""
        raise IOError, error_msg
      end
    end


    def add_sequence( sequence_name, in_msa, add_to_msa )
      seqs = in_msa.find_by_name_start( sequence_name, true )
      if ( seqs.length < 1 )
        error_msg = "sequence \"" + sequence_name + "\" not found in sequence file"
        raise StandardError, error_msg
      end
      if ( seqs.length > 1 )
        error_msg = "sequence \"" + sequence_name + "\" not unique in sequence file"
        raise StandardError, error_msg
      end
      seq = in_msa.get_sequence( seqs[ 0 ] )
      add_to_msa.add_sequence( seq )
    end

    # raises ArgumentError, StandardError
    def extract_domain( sequence,
        number,
        out_of,
        seq_from,
        seq_to,
        in_msa,
        out_msa,
        add_position,
        add_domain_number,
        add_domain_number_as_digit,
        add_domain_number_as_letter,
        trim_name,
        add_species )
      if  number.is_a?( Fixnum ) && ( number < 1 || out_of < 1 || number > out_of )
        error_msg = "impossible: number=" + number.to_s + ", out of=" + out_of.to_s
        raise ArgumentError, error_msg
      end
      if  seq_from < 1 || seq_to < 1 || seq_from >= seq_to
        error_msg = "impossible: seq-f=" + seq_from.to_s + ", seq-t=" + seq_to.to_s
        raise ArgumentError, error_msg
      end
      seqs = in_msa.find_by_name_start( sequence, true )
      if seqs.length < 1
        error_msg = "sequence \"" + sequence + "\" not found in sequence file"
        raise StandardError, error_msg
      end
      if seqs.length > 1
        error_msg = "sequence \"" + sequence + "\" not unique in sequence file"
        raise StandardError, error_msg
      end
      # hmmsearch is 1 based, wheres sequences are 0 bases in this package.
      seq = in_msa.get_sequence( seqs[ 0 ] ).get_subsequence( seq_from - 1, seq_to - 1 )

      orig_name = seq.get_name

      seq.set_name( orig_name.split[ 0 ] )

      if add_position
        seq.set_name( seq.get_name + "_" + seq_from.to_s + "-" + seq_to.to_s )
      end

      if trim_name
        seq.set_name( seq.get_name[ 0, seq.get_name.length - TRIM_BY ] )
      end

      if out_of != 1
        if add_domain_number_as_digit
          seq.set_name( seq.get_name + number.to_s )
        elsif add_domain_number_as_letter
          if number > 25
            error_msg = 'too many identical domains per sequence, cannot use letters to distinguish them'
            raise StandardError, error_msg
          end
          seq.set_name( seq.get_name + ( number + 96 ).chr )
        elsif add_domain_number
          seq.set_name( seq.get_name + "~" + number.to_s + "-" + out_of.to_s )
        end
      end

      # if ( seq.get_name.length > 10 )
      #   error_msg = "sequence name [" + seq.get_name + "] is longer than 10 characters"
      #   raise StandardError, error_msg
      # end

      if add_species
        a = orig_name.rindex "["
        b = orig_name.rindex "]"
        unless a && b
          error_msg = "species not found in " + orig_name
          raise StandardError, error_msg
        end
        species = orig_name[ a .. b ]
        seq.set_name( seq.get_name + " " + species )
      end
      out_msa.add_sequence( seq )
    end

    def is_ignorable?( line )
      return ( line !~ /[A-Za-z0-9-]/ || line =~/^#/ )
    end

  end # class HmmscanDomainExtractor

end # module Evoruby


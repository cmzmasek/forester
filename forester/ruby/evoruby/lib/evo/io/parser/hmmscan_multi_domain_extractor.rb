#
# = lib/evo/io/parser/hmmscan_domain_extractor.rb - HmmscanMultiDomainExtractor class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)
#
# Last modified: 2017/02/20

require 'lib/evo/util/constants'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/parser/hmmscan_parser'

module Evoruby
  class HmmscanMultiDomainExtractor
    def initialize
      @passed = Hash.new
      @non_passsing_domain_architectures = Hash.new
    end

    # raises ArgumentError, IOError, StandardError
    def parse( domain_id,
      hmmscan_output,
      fasta_sequence_file,
      outfile,
      passed_seqs_outfile,
      failed_seqs_outfile,
      e_value_threshold,
      length_threshold,
      add_position,
      add_domain_number,
      add_species,
      log )

      Util.check_file_for_readability( hmmscan_output )
      Util.check_file_for_readability( fasta_sequence_file )
      Util.check_file_for_writability( outfile + ".fasta" )
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
      passed_multi_seqs = Msa.new

      ld = Constants::LINE_DELIMITER

      domain_pass_counter                = 0
      domain_fail_counter                = 0
      passing_domains_per_protein        = 0
      proteins_with_failing_domains      = 0
      domain_not_present_counter         = 0
      protein_counter                    = 1
      max_domain_copy_number_per_protein = -1
      max_domain_copy_number_sequence    = ""
      passing_target_length_sum          = 0
      overall_target_length_sum          = 0
      overall_target_length_min          = 10000000
      overall_target_length_max          = -1
      passing_target_length_min          = 10000000
      passing_target_length_max          = -1

      overall_target_ie_min          = 10000000
      overall_target_ie_max          = -1
      passing_target_ie_min          = 10000000
      passing_target_ie_max          = -1

      hmmscan_datas = []

      hmmscan_parser = HmmscanParser.new( hmmscan_output )
      results = hmmscan_parser.parse

      ####
      # Import: if multiple copies of same domain, threshold need to be same!
      target_domain_ary = Array.new
      target_domain_ary.push(TargetDomain.new('Hexokinase_1', 0.1, -1, 0.5, 1 ))
      target_domain_ary.push(TargetDomain.new('Hexokinase_2', 0.1, -1, 0.5, 1 ))
      target_domain_ary.push(TargetDomain.new('Hexokinase_2', 0.1, -1, 0.5, 1 ))
      target_domain_ary.push(TargetDomain.new('Hexokinase_1', 0.1, -1, 0.5, 1 ))

      #target_domain_ary.push(TargetDomain.new('BH4', 0.1, -1, 0.5, 0 ))
      #target_domain_ary.push(TargetDomain.new('Bcl-2', 0.1, -1, 0.5, 1 ))
      # target_domain_ary.push(TargetDomain.new('Bcl-2_3', 0.01, -1, 0.5, 2 ))

      #  target_domain_ary.push(TargetDomain.new('Nitrate_red_del', 1000, -1, 0.1, 0 ))
      #  target_domain_ary.push(TargetDomain.new('Nitrate_red_del', 1000, -1, 0.1, 1 ))

      #target_domain_ary.push(TargetDomain.new('Chordopox_A33R', 1000, -1, 0.1 ))

      target_domains = Hash.new

      target_domain_architecure = ""

      target_domain_ary.each do |target_domain|
        target_domains[target_domain.name] = target_domain
        if target_domain_architecure.length > 0
          target_domain_architecure += ">"
        end
        target_domain_architecure += target_domain.name
      end

      target_domain_architecure.freeze

      puts target_domain_architecure

      target_domain_names = Set.new

      target_domains.each_key {|key| target_domain_names.add( key ) }

      prev_query_seq_name = nil
      domains_in_query_seq = Array.new
      passing_sequences = Array.new
      total_sequences = 0
      @passed = Hash.new
      results.each do |hmmscan_result|
        if ( prev_query_seq_name != nil ) && ( hmmscan_result.query != prev_query_seq_name )
          if checkit2(domains_in_query_seq, target_domain_names, target_domains, target_domain_architecure)
            passing_sequences.push(domains_in_query_seq)
          end
          domains_in_query_seq = Array.new
          total_sequences += 1
        end
        prev_query_seq_name = hmmscan_result.query
        domains_in_query_seq.push(hmmscan_result)
      end # each result

      if prev_query_seq_name != nil
        total_sequences += 1
        if checkit2(domains_in_query_seq, target_domain_names, target_domains, target_domain_architecure)
          passing_sequences.push(domains_in_query_seq)
        end
      end

      passing_sequences.each do | domains |

        seq = domains[0].query
        # puts seq + "::"
        ################
        if passed_seqs.find_by_name_start( seq, true ).length < 1
          add_sequence( seq, in_msa, passed_multi_seqs )
        else
          error_msg = "this should not have happened"
          raise StandardError, error_msg
        end
        ################

        domains.each do | domain |
          #puts domain.query + ": " + domain.model
        end
        #puts
      end

      puts @non_passsing_domain_architectures

      puts @passed
      puts "total sequences  : " + total_sequences.to_s
      puts "passing sequences: " + passing_sequences.size.to_s

      ####

      prev_query = nil
      saw_target = false

      results.each do | r |

        if ( prev_query != nil ) && ( r.query != prev_query )
          protein_counter += 1
          passing_domains_per_protein = 0
          if !saw_target
            log << domain_not_present_counter.to_s + ": " + prev_query.to_s + " lacks target domain" + ld
            domain_not_present_counter += 1
          end
          saw_target = false
        end

        prev_query = r.query

        if domain_id != r.model
          next
        end

        saw_target = true

        sequence  = r.query
        number    = r.number
        out_of    = r.out_of
        env_from  = r.env_from
        env_to    = r.env_to
        i_e_value = r.i_e_value
        prev_query = r.query

        length = 1 + env_to - env_from

        overall_target_length_sum += length
        if length > overall_target_length_max
          overall_target_length_max = length
        end
        if length < overall_target_length_min
          overall_target_length_min = length
        end

        if i_e_value > overall_target_ie_max
          overall_target_ie_max = i_e_value
        end
        if i_e_value < overall_target_ie_min
          overall_target_ie_min = i_e_value
        end

        if ( ( ( e_value_threshold < 0.0 ) || ( i_e_value <= e_value_threshold ) ) &&
        ( ( length_threshold <= 0 ) || ( length >= length_threshold.to_f ) ) )
          hmmscan_datas << HmmscanData.new( sequence, number, out_of, env_from, env_to, i_e_value )
          passing_target_length_sum += length
          passing_domains_per_protein += 1
          if length > passing_target_length_max
            passing_target_length_max = length
          end
          if length < passing_target_length_min
            passing_target_length_min = length
          end
          if i_e_value > passing_target_ie_max
            passing_target_ie_max = i_e_value
          end
          if i_e_value < passing_target_ie_min
            passing_target_ie_min = i_e_value
          end
          if ( passing_domains_per_protein > max_domain_copy_number_per_protein )
            max_domain_copy_number_sequence    = sequence
            max_domain_copy_number_per_protein = passing_domains_per_protein
          end
        else # no pass
          log << domain_fail_counter.to_s + ": " + sequence.to_s + " fails threshold(s)"
          if ( ( e_value_threshold.to_f >= 0.0 ) && ( i_e_value > e_value_threshold ) )
            log << " iE=" + i_e_value.to_s
          end
          if ( ( length_threshold.to_f > 0 ) && ( env_to - env_from + 1 ) < length_threshold.to_f )
            le = env_to - env_from + 1
            log << " l=" + le.to_s
          end
          log << ld
          domain_fail_counter += 1
        end

        if number > out_of
          error_msg = "number > out_of (this should not have happened)"
          raise StandardError, error_msg
        end

        if number == out_of
          if !hmmscan_datas.empty?
            process_hmmscan_datas( hmmscan_datas,
            in_msa,
            add_domain_number,
            out_msa )
            domain_pass_counter += hmmscan_datas.length
            if passed_seqs.find_by_name_start( sequence, true ).length < 1
              add_sequence( sequence, in_msa, passed_seqs )
            else
              error_msg = "this should not have happened"
              raise StandardError, error_msg
            end
          else # no pass
            if failed_seqs.find_by_name_start( sequence, true ).length < 1
              add_sequence( sequence, in_msa, failed_seqs )
              proteins_with_failing_domains += 1
            else
              error_msg = "this should not have happened"
              raise StandardError, error_msg
            end
          end
          hmmscan_datas.clear
        end

      end # results.each do | r |

      if (prev_query != nil) && (!saw_target)
        log << domain_not_present_counter.to_s + ": " + prev_query.to_s + " lacks target domain" + ld
        domain_not_present_counter += 1
      end

      if domain_pass_counter < 1
        error_msg = "no domain sequences were extracted"
        raise IOError, error_msg
      end

      if ( domain_not_present_counter + passed_seqs.get_number_of_seqs + proteins_with_failing_domains ) != protein_counter
        error_msg = "not present + passing + not passing != proteins in sequence (fasta) file (this should not have happened)"
        raise StandardError, error_msg
      end

      puts
      log << ld

      log << ld
      avg_pass = ( passing_target_length_sum / domain_pass_counter )
      puts( "Passing target domain lengths: average: " + avg_pass.to_s  )
      log << "Passing target domain lengths: average: " + avg_pass.to_s
      log << ld
      puts( "Passing target domain lengths: min-max: " + passing_target_length_min.to_s + " - "  + passing_target_length_max.to_s)
      log << "Passing target domain lengths: min-max: " + passing_target_length_min.to_s + " - "  + passing_target_length_max.to_s
      log << ld
      puts( "Passing target domain iE:      min-max: " + passing_target_ie_min.to_s + " - "  + passing_target_ie_max.to_s)
      log << "Passing target domain iE:      min-max: " + passing_target_ie_min.to_s + " - "  + passing_target_ie_max.to_s
      log << ld
      puts( "Passing target domains:            sum: " + domain_pass_counter.to_s  )
      log << "Passing target domains:            sum: " + domain_pass_counter.to_s
      log << ld
      log << ld
      puts
      sum = domain_pass_counter + domain_fail_counter
      avg_all = overall_target_length_sum / sum
      puts( "All target domain lengths:     average: " + avg_all.to_s  )
      log << "All target domain lengths:     average: " + avg_all.to_s
      log << ld
      puts( "All target domain lengths:     min-max: " + overall_target_length_min.to_s + " - "  + overall_target_length_max.to_s)
      log << "All target domain lengths:     min-max: " + overall_target_length_min.to_s + " - "  + overall_target_length_max.to_s
      log << ld
      puts( "All target target domain iE:   min-max: " + overall_target_ie_min.to_s + " - "  + overall_target_ie_max.to_s)
      log << "All target target domain iE:   min-max: " + overall_target_ie_min.to_s + " - "  + overall_target_ie_max.to_s
      log << ld
      puts( "All target domains:                sum: " + sum.to_s  )
      log << "All target domains:                sum: " + sum.to_s

      puts
      puts( "Proteins with passing target domain(s): " + passed_multi_seqs.get_number_of_seqs.to_s )
      #puts( "Proteins with passing target domain(s): " + passed_seqs.get_number_of_seqs.to_s )
      puts( "Proteins with no passing target domain: " + proteins_with_failing_domains.to_s )
      puts( "Proteins with no target domain        : " + domain_not_present_counter.to_s )

      log << ld
      log << ld
      puts
      puts( "Max target domain copy number per protein: " + max_domain_copy_number_per_protein.to_s )
      log << "Max target domain copy number per protein: " + max_domain_copy_number_per_protein.to_s
      log << ld

      if ( max_domain_copy_number_per_protein > 1 )
        puts( "First target protein with this copy number: " + max_domain_copy_number_sequence )
        log << "First target protein with this copy number: " + max_domain_copy_number_sequence
        log << ld
      end

      write_msa( out_msa, outfile + ".fasta"  )
      write_msa( passed_multi_seqs, passed_seqs_outfile )
      write_msa( failed_seqs, failed_seqs_outfile )

      log << ld
      log << "passing target domains                       : " + domain_pass_counter.to_s + ld
      log << "failing target domains                       : " + domain_fail_counter.to_s + ld
      log << "proteins in sequence (fasta) file            : " + in_msa.get_number_of_seqs.to_s + ld
      log << "proteins in hmmscan outputfile               : " + protein_counter.to_s + ld
      log << "proteins with passing target domain(s)       : " + passed_seqs.get_number_of_seqs.to_s + ld
      log << "proteins with no passing target domain       : " + proteins_with_failing_domains.to_s + ld
      log << "proteins with no target domain               : " + domain_not_present_counter.to_s + ld

      log << ld

      return domain_pass_counter

    end # parse

    private

    def checkit2(domains_in_query_seq, target_domain_names, target_domains, target_domain_architecture = nil)
      domain_names_in_query_seq = Set.new
      domains_in_query_seq.each do |domain|
        domain_names_in_query_seq.add(domain.model)
      end
      if (target_domain_names.subset?(domain_names_in_query_seq))

        passed_domains = Array.new

        domains_in_query_seq.each do |domain|
          if target_domains.has_key?(domain.model)
            target_domain = target_domains[domain.model]
            if (target_domain.i_e_value != nil) && (target_domain.i_e_value >= 0)
              if domain.i_e_value > target_domain.i_e_value
                return false
              end
            end
            if (target_domain.abs_len != nil) && (target_domain.abs_len > 0)
              length = 1 + domain.env_to - domain.env_from
              if length < target_domain.abs_len
                return false
              end
            end
            if (target_domain.rel_len != nil) && (target_domain.rel_len > 0)
              length = 1 + domain.env_to - domain.env_from
              if length < (target_domain.rel_len * domain.tlen)
                return false
              end
            end

            # For domain architecture comparison
            if target_domain_architecture != nil
              passed_domains.push(domain)
            end

            if (@passed.key?(domain.model))
              @passed[domain.model] =  @passed[domain.model] + 1
            else
              @passed[domain.model] = 1
            end
          end
        end

      else
        return false
      end
      # Compare architectures
      if target_domain_architecture != nil
        if !compareArchitectures(target_domain_architecture, passed_domains)
          return false
        end
      end
      true
    end

    def compareArchitectures(target_domain_architecture, passed_domains)
      passed_domains.sort! { |a,b| a.ali_from <=> b.ali_from }
      domain_architecture = ""
      passed_domains.each do |domain|
        if domain_architecture.length > 0
          domain_architecture += ">"
        end
        domain_architecture += domain.model
      end
      pass = (target_domain_architecture == domain_architecture)
      if !pass
        if (@non_passsing_domain_architectures.key?(domain_architecture))
          @non_passsing_domain_architectures[domain_architecture] = @non_passsing_domain_architectures[domain_architecture] + 1
        else
          @non_passsing_domain_architectures[domain_architecture] = 1
        end
      end
      return pass
    end

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

    def process_hmmscan_datas( hmmscan_datas,
      in_msa,
      add_domain_number,
      out_msa )

      actual_out_of = hmmscan_datas.size

      seq_name = ""
      prev_seq_name = nil

      hmmscan_datas.each_with_index do |hmmscan_data, index|
        if hmmscan_data.number < ( index + 1 )
          error_msg = "hmmscan_data.number < ( index + 1 ) (this should not have happened)"
          raise StandardError, error_msg
        end

        seq_name = hmmscan_data.seq_name

        extract_domain( seq_name,
        index + 1,
        actual_out_of,
        hmmscan_data.env_from,
        hmmscan_data.env_to,
        in_msa,
        out_msa,
        add_domain_number )

        if prev_seq_name && prev_seq_name != seq_name
          error_msg = "this should not have happened"
          raise StandardError, error_msg
        end
        prev_seq_name = seq_name
      end # each
    end # def process_hmmscan_data

    def extract_domain( sequence,
      number,
      out_of,
      seq_from,
      seq_to,
      in_msa,
      out_msa,
      add_domain_number)
      if number.is_a?( Fixnum ) && ( number < 1 || out_of < 1 || number > out_of )
        error_msg = "number=" + number.to_s + ", out of=" + out_of.to_s
        raise StandardError, error_msg
      end
      if seq_from < 1 || seq_to < 1 || seq_from >= seq_to
        error_msg = "impossible: seq-from=" + seq_from.to_s + ", seq-to=" + seq_to.to_s
        raise StandardError, error_msg
      end
      seqs = in_msa.find_by_name_start( sequence, true )
      if seqs.length < 1
        error_msg = "sequence \"" + sequence + "\" not found in sequence file"
        raise IOError, error_msg
      end
      if seqs.length > 1
        error_msg = "sequence \"" + sequence + "\" not unique in sequence file"
        raise IOError, error_msg
      end
      # hmmscan is 1 based, whereas sequences are 0 bases in this package.
      seq = in_msa.get_sequence( seqs[ 0 ] ).get_subsequence( seq_from - 1, seq_to - 1 )

      orig_name = seq.get_name

      seq.set_name( orig_name.split[ 0 ] )

      if out_of != 1 && add_domain_number
        seq.set_name( seq.get_name + "~" + number.to_s + "-" + out_of.to_s )
      end

      out_msa.add_sequence( seq )
    end

    def is_ignorable?( line )
      return ( line !~ /[A-Za-z0-9-]/ || line =~/^#/ )
    end

  end # class HmmscanDomainExtractor

  class HmmscanData
    def initialize( seq_name, number, out_of, env_from, env_to, i_e_value )
      @seq_name = seq_name
      @number = number
      @out_of = out_of
      @env_from = env_from
      @env_to = env_to
      @i_e_value = i_e_value
    end
    attr_reader :seq_name, :number, :out_of, :env_from, :env_to, :i_e_value
  end

  class TargetDomain
    def initialize( name, i_e_value, abs_len, rel_len, position )
      @name = name
      @i_e_value = i_e_value
      @abs_len = abs_len
      @percent_len = rel_len
      @position = position
    end
    attr_reader :name, :i_e_value, :abs_len, :rel_len, :position
  end

end # module Evoruby


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

      hmmscan_parser = HmmscanParser.new( hmmscan_output )
      results = hmmscan_parser.parse

      ####
      # Import: if multiple copies of same domain, threshold need to be same!
      target_domain_ary = Array.new
      target_domain_ary.push(TargetDomain.new('Hexokinase_1', 1e-6, -1, 0.6, 1 ))
      target_domain_ary.push(TargetDomain.new('Hexokinase_2', 1e-6, -1, 0.6, 1 ))
      # target_domain_ary.push(TargetDomain.new('Hexokinase_2', 0.1, -1, 0.5, 1 ))
      # target_domain_ary.push(TargetDomain.new('Hexokinase_1', 0.1, -1, 0.5, 1 ))

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

      puts 'Target domain architecture: ' + target_domain_architecure

      target_domain_names = Set.new

      target_domains.each_key {|key| target_domain_names.add( key ) }

      prev_query_seq_name = nil
      domains_in_query_seq = Array.new
      passing_sequences = Array.new
      total_sequences = 0
      @passed = Hash.new
      out_domain_msas = Hash.new
      out_domain_architecture_msa = Msa.new
      results.each do |hmmscan_result|
        if ( prev_query_seq_name != nil ) && ( hmmscan_result.query != prev_query_seq_name )
          if checkit2(domains_in_query_seq, target_domain_names, target_domains, in_msa, out_domain_msas, out_domain_architecture_msa, target_domain_architecure)
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
        if checkit2(domains_in_query_seq, target_domain_names, target_domains, in_msa, out_domain_msas, out_domain_architecture_msa, target_domain_architecure)
          passing_sequences.push(domains_in_query_seq)
        end
      end

      out_domain_msas.keys.sort.each do |domain_name|
        puts "writing #{domain_name}"
        write_msa( out_domain_msas[domain_name], domain_name + ".fasta" )
      end

      puts "writing target_domain_architecure"
      write_msa( out_domain_architecture_msa, "target_domain_architecure" + ".fasta" )

      passing_sequences.each do | domains |

        seq = domains[0].query
        # puts seq + "::"

        if passed_seqs.find_by_name_start( seq, true ).length < 1
          add_sequence( seq, in_msa, passed_multi_seqs )
        else
          error_msg = "this should not have happened"
          raise StandardError, error_msg
        end

        domains.each do | domain |
          #puts domain.query + ": " + domain.model
        end
        #puts
      end

      puts

      puts 'Non passing domain architectures:'
      @non_passsing_domain_architectures = @non_passsing_domain_architectures.sort{|a, b|a<=>b}.to_h
      @non_passsing_domain_architectures.each do |da, count|
        puts da + ': ' + count.to_s
      end

      puts

      puts 'Passing domain counts:'
      @passed = @passed.sort{|a, b|a<=>b}.to_h
      @passed.each do |dom, count|
        puts dom + ': ' + count.to_s
      end

      puts

      puts "total sequences  : " + total_sequences.to_s
      puts "passing sequences: " + passing_sequences.size.to_s

      # write_msa( out_msa, outfile + ".fasta"  )
      # write_msa( passed_multi_seqs, passed_seqs_outfile )
      # write_msa( failed_seqs, failed_seqs_outfile )

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

    # domains_in_query_seq: Array of HmmscanResult
    # target_domain_names: Set of String
    # target_domains: Hash String->TargetDomain
    # target_domain_architecture: String
    def checkit2(domains_in_query_seq,
      target_domain_names,
      target_domains,
      in_msa,
      out_domain_msas,
      out_domain_architecture_msa,
      target_domain_architecture)

      domain_names_in_query_seq = Set.new
      domains_in_query_seq.each do |domain|
        domain_names_in_query_seq.add(domain.model)
      end
      if (target_domain_names.subset?(domain_names_in_query_seq))

        passed_domains = Array.new
        passed_domains_counts = Hash.new

        domains_in_query_seq.each do |domain|
          if target_domains.has_key?(domain.model)
            target_domain = target_domains[domain.model]

            if (target_domain.i_e_value != nil) && (target_domain.i_e_value >= 0)
              if domain.i_e_value > target_domain.i_e_value
                next
              end
            end
            if (target_domain.abs_len != nil) && (target_domain.abs_len > 0)
              length = 1 + domain.env_to - domain.env_from
              if length < target_domain.abs_len
                next
              end
            end
            if (target_domain.rel_len != nil) && (target_domain.rel_len > 0)
              length = 1 + domain.env_to - domain.env_from
              puts (target_domain.rel_len * domain.tlen)
              if length < (target_domain.rel_len * domain.tlen)
                next
              end
            end

            passed_domains.push(domain)

            if (passed_domains_counts.key?(domain.model))
              passed_domains_counts[domain.model] = passed_domains_counts[domain.model] + 1
            else
              passed_domains_counts[domain.model] = 1
            end

            if (@passed.key?(domain.model))
              @passed[domain.model] = @passed[domain.model] + 1
            else
              @passed[domain.model] = 1
            end
          end # if target_domains.has_key?(domain.model)
        end # domains_in_query_seq.each do |domain|
      else
        return false
      end

      passed_domains.sort! { |a,b| a.ali_from <=> b.ali_from }
      # Compare architectures
      if !compareArchitectures(target_domain_architecture, passed_domains)
        return false
      end

      domain_counter = Hash.new

      min_env_from = 10000000
      max_env_from = 0
      min_env_to = 10000000
      max_env_to = 0

      query_seq = nil

      concatenated_domains = ''

      passed_domains.each do |domain|
        domain_name = domain.model

        unless out_domain_msas.has_key? domain_name
          out_domain_msas[ domain_name ] = Msa.new
        end

        if domain_counter.key?(domain_name)
          domain_counter[domain_name] = domain_counter[domain_name] + 1
        else
          domain_counter[domain_name] = 1
        end

        if query_seq == nil
          query_seq = domain.query
          query_seq.freeze
        elsif query_seq != domain.query
          error_msg = "seq names do not match: this should not have happened"
          raise StandardError, error_msg
        end
        puts "query: " + query_seq

        extracted_domain_seq = extract_domain( query_seq,
        domain_counter[domain_name],
        passed_domains_counts[domain_name],
        domain.env_from,
        domain.env_to,
        in_msa,
        out_domain_msas[ domain_name ] )

        concatenated_domains += extracted_domain_seq

        if domain.env_from < min_env_from
          min_env_from = domain.env_from
        end
        if domain.env_from > max_env_from
          max_env_from = domain.env_from
        end
        if domain.env_to < min_env_to
          min_env_to = domain.env_to
        end
        if domain.env_to > max_env_to
          max_env_to = domain.env_to
        end
      end

      puts "env from: #{min_env_from} - #{max_env_from}"
      puts "env to  : #{min_env_to} - #{max_env_to}"

      extract_domain( query_seq,
      1,
      1,
      min_env_from,
      max_env_to,
      in_msa,
      out_domain_architecture_msa )

      puts passed_domains_counts
      puts concatenated_domains
      true
    end

    # passed_domains needs to be sorted!
    def compareArchitectures(target_domain_architecture, passed_domains, strict = false)
      domain_architecture = ""
      passed_domains.each do |domain|
        if domain_architecture.length > 0
          domain_architecture += ">"
        end
        domain_architecture += domain.model
      end
      pass = false
      if strict
        pass = (target_domain_architecture == domain_architecture)
      else
        pass = (domain_architecture.index(target_domain_architecture) != nil)
      end
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
      File.delete(filename) if File.exist?(filename) #TODO remove me
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

    def extract_domain( seq_name,
      number,
      out_of,
      seq_from,
      seq_to,
      in_msa,
      out_msa)
      if number < 1 || out_of < 1 || number > out_of
        error_msg = "number=" + number.to_s + ", out of=" + out_of.to_s
        raise StandardError, error_msg
      end
      if seq_from < 1 || seq_to < 1 || seq_from >= seq_to
        error_msg = "impossible: seq-from=" + seq_from.to_s + ", seq-to=" + seq_to.to_s
        raise StandardError, error_msg
      end
      seqs = in_msa.find_by_name_start(seq_name, true)
      if seqs.length < 1
        error_msg = "sequence \"" + seq_name + "\" not found in sequence file"
        raise IOError, error_msg
      end
      if seqs.length > 1
        error_msg = "sequence \"" + seq_name + "\" not unique in sequence file"
        raise IOError, error_msg
      end
      # hmmscan is 1 based, whereas sequences are 0 bases in this package.
      seq = in_msa.get_sequence( seqs[ 0 ] ).get_subsequence( seq_from - 1, seq_to - 1 )

      orig_name = seq.get_name

      seq.set_name( orig_name.split[ 0 ] )

      if out_of != 1
        seq.set_name( seq.get_name + "~" + number.to_s + "-" + out_of.to_s )
      end

      out_msa.add_sequence( seq )

      seq.get_sequence_as_string
    end

  end # class HmmscanMultiDomainExtractor

  class TargetDomain
    def initialize( name, i_e_value, abs_len, rel_len, position )
      if (name == nil) || name.size < 1
        error_msg = "target domain name nil or empty"
        raise StandardError, error_msg
      end
      if rel_len > 1
        error_msg = "target domain relative length is greater than 1"
        raise StandardError, error_msg
      end
      @name = name
      @i_e_value = i_e_value
      @abs_len = abs_len
      @rel_len = rel_len
      @position = position
    end
    attr_reader :name, :i_e_value, :abs_len, :rel_len, :position
  end

end # module Evoruby


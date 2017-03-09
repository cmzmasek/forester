#
# = lib/evo/io/parser/hmmscan_domain_extractor.rb - HmmscanMultiDomainExtractor class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)
#
# Last modified: 2017/03/08

####
# Import: if multiple copies of same domain, thresholds need to be same!
####

require 'lib/evo/util/constants'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/parser/hmmscan_parser'

module Evoruby
  class HmmscanMultiDomainExtractor

    TESTING = false
    OUTPUT_ID = 'mdsx'
    DOMAIN_DELIMITER = ' -- '

    IE_CUTOFF_FOR_DA_OVERVIEW      = 1E-6
    REL_LEN_CUTOFF_FOR_DA_OVERVIEW = 0.5

    PASSING_FL_SEQS_SUFFIX    = "_#{OUTPUT_ID}_passing_full_length_seqs.fasta"
    FAILING_FL_SEQS_SUFFIX    = "_#{OUTPUT_ID}_failing_full_length_seqs.fasta"
    TARGET_DA_SUFFIX          = "_#{OUTPUT_ID}_target_da.fasta"
    CONCAT_TARGET_DOM_SUFFIX  = "_#{OUTPUT_ID}_concat_target_doms.fasta"
    TARGET_DOM_OUTPUT_MIDPART = "_#{OUTPUT_ID}_target_dom_"
    LOG_FILE_SUFFIX           = "_#{OUTPUT_ID}.log"
    def initialize
      @passing_domains_data = nil
      @failing_domains_data = nil
      @failing_domain_architectures = nil
      @passsing_domain_architectures = nil
      @failing_proteins_bc_not_all_target_doms_present = nil
      @failing_proteins_bc_missing_cutoffs = nil
      @encountered_domain_architectures = nil
      @log = nil
    end

    # raises ArgumentError, IOError, StandardError
    def parse( target_da,
      hmmscan_output,
      fasta_sequence_file,
      outfile_base,
      log_str )

      passing_fl_seqs_outfile   = outfile_base + PASSING_FL_SEQS_SUFFIX
      failing_fl_seqs_outfile   = outfile_base + FAILING_FL_SEQS_SUFFIX
      target_da_outfile         = outfile_base + TARGET_DA_SUFFIX
      concat_target_dom_outfile = outfile_base + CONCAT_TARGET_DOM_SUFFIX
      logfile                   = outfile_base + LOG_FILE_SUFFIX

      Util.check_file_for_readability( hmmscan_output )
      Util.check_file_for_readability( fasta_sequence_file )
      Util.check_file_for_writability( passing_fl_seqs_outfile )
      Util.check_file_for_writability( failing_fl_seqs_outfile )
      Util.check_file_for_writability( target_da_outfile )
      Util.check_file_for_writability( concat_target_dom_outfile )
      Util.check_file_for_writability( logfile )

      @log = log_str

      in_msa = nil
      factory = MsaFactory.new()
      in_msa = factory.create_msa_from_file( fasta_sequence_file, FastaParser.new() )

      if ( in_msa == nil || in_msa.get_number_of_seqs() < 1 )
        error_msg = "could not find fasta sequences in " + fasta_sequence_file
        raise IOError, error_msg
      end

      failed_seqs_msa = Msa.new
      passed_seqs_msa = Msa.new

      hmmscan_parser = HmmscanParser.new( hmmscan_output )
      results = hmmscan_parser.parse

      target_domain_ary = parse_da target_da

      target_domains = Hash.new

      target_domain_architecure = ''

      target_domain_ary.each do |target_domain|
        target_domains[target_domain.name] = target_domain
        if target_domain_architecure.length > 0
          target_domain_architecure << DOMAIN_DELIMITER
        end
        target_domain_architecure << target_domain.name
      end

      target_domain_architecure.freeze

      log 'Hmmscan outputfile             : ' + hmmscan_output
      log 'Full length fasta sequence file: ' + fasta_sequence_file
      log 'Target domain architecture     : ' + target_domain_architecure
      target_domain_ary.each do |x|
        log x.to_str
      end

      target_domain_names = Set.new

      target_domains.each_key {|key| target_domain_names.add( key ) }

      prev_query_seq_name = nil
      domains_in_query_seq = Array.new
      passing_sequences = Array.new # This will be an Array of Array of HmmscanResult for the same query seq.
      failing_sequences = Array.new # This will be an Array of Array of HmmscanResult for the same query seq.
      total_sequences = 0

      @failing_domain_architectures = Hash.new
      @passsing_domain_architectures = Hash.new
      @passing_domains_data = Hash.new
      @failing_domains_data = Hash.new
      @encountered_domain_architectures= Hash.new
      @failing_proteins_bc_not_all_target_doms_present = 0
      @failing_proteins_bc_missing_cutoffs = 0
      out_domain_msas = Hash.new
      out_domain_architecture_msa = Msa.new
      out_concatenated_domains_msa = Msa.new
      results.each do |hmmscan_result|
        if ( prev_query_seq_name != nil ) && ( hmmscan_result.query != prev_query_seq_name )
          if compare(domains_in_query_seq, target_domain_names, target_domains, in_msa, out_domain_msas, out_domain_architecture_msa, out_concatenated_domains_msa, target_domain_architecure)
            passing_sequences.push(domains_in_query_seq)
          else
            failing_sequences.push(domains_in_query_seq)
          end
          domains_in_query_seq = Array.new
          total_sequences += 1
        end
        prev_query_seq_name = hmmscan_result.query
        domains_in_query_seq.push(hmmscan_result)
      end # each result

      if prev_query_seq_name != nil
        if compare(domains_in_query_seq, target_domain_names, target_domains, in_msa, out_domain_msas, out_domain_architecture_msa, out_concatenated_domains_msa, target_domain_architecure)
          passing_sequences.push(domains_in_query_seq)
        else
          failing_sequences.push(domains_in_query_seq)
        end
        total_sequences += 1
      end

      if in_msa.get_number_of_seqs < total_sequences
        error_msg = "hmmscan output contains more protein sequences than fasta sequence file"
        raise IOError, error_msg
      end

      log
      log 'Domain architecture overview using default iE-cutoff ' +
      IE_CUTOFF_FOR_DA_OVERVIEW.to_s + ' and relative length cutoff ' + REL_LEN_CUTOFF_FOR_DA_OVERVIEW.to_s + ':'

      @encountered_domain_architectures = @encountered_domain_architectures.sort_by {|k, v| v}.reverse.to_h
      counter = 1;
      @encountered_domain_architectures.each do |k, v|
        log counter.to_s.rjust(2) + ') ' +  v.to_s.rjust(5) + ': ' + k
        counter += 1
        if counter > 40
          break
        end
      end

      log
      log 'Passing domain arrangements of target domain(s):'
      @passsing_domain_architectures = @passsing_domain_architectures.sort{|a, b|a<=>b}.to_h
      passing_da_sum = 0
      @passsing_domain_architectures.each do |da, count|
        passing_da_sum += count
        log count.to_s.rjust(4) + ': ' + da
      end
      log
      log 'Failing domain arrangements of target domain(s):'
      @failing_domain_architectures = @failing_domain_architectures.sort{|a, b|a<=>b}.to_h
      failing_da_sum = 0
      @failing_domain_architectures .each do |da, count|
        failing_da_sum += count
        log count.to_s.rjust(4) + ': ' + da
      end
      log
      log 'Passing target domain(s):'
      @passing_domains_data = @passing_domains_data.sort{|a, b|a<=>b}.to_h
      @passing_domains_data.each do |n, d|
        log d.to_str
      end
      log
      log'Failing target domain(s):'
      @failing_domains_data = @failing_domains_data.sort{|a, b|a<=>b}.to_h
      @failing_domains_data.each do |n, d|
        log d.to_str
      end

      unless total_sequences == (passing_sequences.size + failing_sequences.size)
        error_msg = "this should not have happened: total seqs not equal to passing plus failing seqs"
        raise StandardError, error_msg
      end

      unless failing_sequences.size == (@failing_proteins_bc_not_all_target_doms_present + @failing_proteins_bc_missing_cutoffs)
        error_msg = "this should not have happened: failing seqs sums not consistent"
        raise StandardError, error_msg
      end

      unless @failing_proteins_bc_missing_cutoffs >= failing_da_sum
        error_msg = "this should not have happened: failing seqs larger than failing da sum"
        raise StandardError, error_msg
      end

      unless passing_sequences.size == passing_da_sum
        error_msg = "this should not have happened: passing seqs not equal to passing da sum"
        raise StandardError, error_msg
      end

      log
      log "Protein sequences in sequence (fasta) file: " + in_msa.get_number_of_seqs.to_s.rjust(5)
      log "Protein sequences in hmmscan output file  : " + total_sequences.to_s.rjust(5)
      log "  Passing protein sequences               : " + passing_sequences.size.to_s.rjust(5)
      log "  Failing protein sequences               : " + failing_sequences.size.to_s.rjust(5)
      log "    Not all target domain present         : " + @failing_proteins_bc_not_all_target_doms_present.to_s.rjust(5)
      log "    Target domain(s) failing cutoffs      : " + @failing_proteins_bc_missing_cutoffs.to_s.rjust(5)
      log

      out_domain_msas.keys.sort.each do |domain_name|
        file_name = outfile_base + TARGET_DOM_OUTPUT_MIDPART + domain_name + '.fasta'
        write_msa out_domain_msas[domain_name], file_name
        log "Wrote passing target domain sequence for " +  domain_name.ljust(16) + ': ' + file_name
      end

      write_msa out_domain_architecture_msa, target_da_outfile
      log 'Wrote target domain architecture                         : ' + target_da_outfile

      write_msa out_concatenated_domains_msa, concat_target_dom_outfile
      log 'Wrote concatenated target domain(s)                      : ' + concat_target_dom_outfile

      passing_sequences.each do | domains |
        query_name = domains[0].query
        if (!TESTING) || (passed_seqs_msa.find_by_name_start( query_name, true ).length < 1)
          add_sequence( query_name, in_msa, passed_seqs_msa )
        else
          error_msg = 'this should not have happened'
          raise StandardError, error_msg
        end
      end

      failing_sequences.each do | domains |
        query_name = domains[0].query
        if (!TESTING) || (failed_seqs_msa.find_by_name_start( query_name, true ).length < 1)
          add_sequence( query_name, in_msa, failed_seqs_msa )
        else
          error_msg = 'this should not have happened'
          raise StandardError, error_msg
        end
      end

      write_msa passed_seqs_msa, passing_fl_seqs_outfile
      log 'Wrote passing full length protein sequences              : ' + passing_fl_seqs_outfile
      write_msa failed_seqs_msa, failing_fl_seqs_outfile
      log 'Wrote failing full length protein sequences              : ' + failing_fl_seqs_outfile

      begin
        f = File.open( logfile, 'w' )
        f.print( @log )
        f.close
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end
      log 'Wrote log file                                           : ' + logfile

    end # parse

    private

    # domains: Array of HmmscanResult
    def collect(domains, ie_cutoff, rel_len_cutoff)
      passed = Array.new
      domains.each do |domain|
        length = 1 + domain.env_to - domain.env_from
        if (domain.i_e_value <= ie_cutoff) && (length >= (rel_len_cutoff * domain.tlen))
          passed.push domain
        end
      end
      passed.sort! { |a,b| a.ali_from <=> b.ali_from }
      s = ''
      passed.each do |domain|
        if s.length > 0
          s << DOMAIN_DELIMITER
        end
        s << domain.model
      end
      if s.length > 0
        if @encountered_domain_architectures.has_key? s
          @encountered_domain_architectures[s] = 1 + @encountered_domain_architectures[s]
        else
          @encountered_domain_architectures[s] = 1
        end
      end
    end

    # domains_in_query_seq: Array of HmmscanResult
    # target_domain_names: Set of String
    # target_domains: Hash String->TargetDomain
    # target_domain_architecture: String
    def compare(domains_in_query_seq,
      target_domain_names,
      target_domains,
      in_msa,
      out_domain_msas,
      out_domain_architecture_msa,
      out_contactenated_domains_msa,
      target_domain_architecture)

      collect(domains_in_query_seq, IE_CUTOFF_FOR_DA_OVERVIEW, REL_LEN_CUTOFF_FOR_DA_OVERVIEW)

      domain_names_in_query_seq = Set.new
      domains_in_query_seq.each do |domain|
        domain_names_in_query_seq.add(domain.model)
      end

      passed_domains = Array.new
      passed_domains_counts = Hash.new

      if (domain_names_in_query_seq.length > 0) && (target_domain_names.subset?(domain_names_in_query_seq))

        domains_in_query_seq.each do |domain|
          if target_domains.has_key?(domain.model)
            target_domain = target_domains[domain.model]

            if target_domain.i_e_value >= 0
              if domain.i_e_value > target_domain.i_e_value
                addToFailingDomainData(domain)
                next
              end
            end
            if target_domain.abs_len > 0
              length = 1 + domain.env_to - domain.env_from
              if length < target_domain.abs_len
                addToFailingDomainData(domain)
                next
              end
            end
            if target_domain.rel_len > 0
              length = 1 + domain.env_to - domain.env_from
              if length < (target_domain.rel_len * domain.tlen)
                addToFailingDomainData(domain)
                next
              end
            end

            passed_domains.push(domain)

            if (passed_domains_counts.key?(domain.model))
              passed_domains_counts[domain.model] = passed_domains_counts[domain.model] + 1
            else
              passed_domains_counts[domain.model] = 1
            end

            addToPassingDomainData(domain)

          end # if target_domains.has_key?(domain.model)
        end # domains_in_query_seq.each do |domain|
      else
        @failing_proteins_bc_not_all_target_doms_present += 1
        return false
      end # target_domain_names.subset?(domain_names_in_query_seq)

      if passed_domains.length < 1
        @failing_proteins_bc_missing_cutoffs += 1
        return false
      end

      passed_domains.sort! { |a,b| a.ali_from <=> b.ali_from }
      # Compare architectures
      if !compareArchitectures(target_domain_architecture, passed_domains, false)
        @failing_proteins_bc_missing_cutoffs += 1
        return false
      end

      domain_counter = Hash.new

      min_env_from = 10000000
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
        if domain.env_to > max_env_to
          max_env_to = domain.env_to
        end
      end

      extract_domain( query_seq,
      1,
      1,
      min_env_from,
      max_env_to,
      in_msa,
      out_domain_architecture_msa )

      out_contactenated_domains_msa.add( query_seq, concatenated_domains)

      true
    end

    def addToPassingDomainData(domain)
      unless ( @passing_domains_data.key?(domain.model))
        @passing_domains_data[domain.model] = DomainData.new(domain.model)
      end
      @passing_domains_data[domain.model].add( domain.model, (1 + domain.env_to - domain.env_from), domain.i_e_value)
    end

    def addToFailingDomainData(domain)
      unless ( @failing_domains_data.key?(domain.model))
        @failing_domains_data[domain.model] = DomainData.new(domain.model)
      end
      @failing_domains_data[domain.model].add( domain.model, (1 + domain.env_to - domain.env_from), domain.i_e_value)
    end

    # passed_domains needs to be sorted!
    def compareArchitectures(target_domain_architecture, passed_domains, strict = false)
      domain_architecture = ''
      passed_domains.each do |domain|
        if domain_architecture.length > 0
          domain_architecture += DOMAIN_DELIMITER
        end
        domain_architecture += domain.model
      end
      pass = false
      if strict
        pass = (target_domain_architecture == domain_architecture)
      else
        pass = (domain_architecture.index(target_domain_architecture) != nil)
      end
      if pass
        if (@passsing_domain_architectures.key?(domain_architecture))
          @passsing_domain_architectures[domain_architecture] = @passsing_domain_architectures[domain_architecture] + 1
        else
          @passsing_domain_architectures[domain_architecture] = 1
        end
      else
        if ( @failing_domain_architectures.key?(domain_architecture))
          @failing_domain_architectures[domain_architecture] =  @failing_domain_architectures[domain_architecture] + 1
        else
          @failing_domain_architectures[domain_architecture] = 1
        end
      end
      return pass
    end

    def write_msa( msa, filename )
      io = MsaIO.new()
      w = FastaWriter.new()
      w.set_line_width( 60 )
      w.clean( true )
      File.delete(filename) if File.exist?(filename)
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

    def parse_da( target_da_str )
      target_domain_hash = Hash.new
      target_domain_ary = Array.new
      target_das = target_da_str.split '--'
      target_das.each do |x|
        inds = x.split '='
        unless inds.size == 4
          raise IOError, 'domain architecture is ill formatted: ' + x
        end
        target_domain_name = inds[0]
        ie_cutoff = Float(inds[1])
        abs_len_cutoff = Integer(inds[2])
        rel_len_cutoff = Float(inds[3])
        if target_domain_hash.has_key? target_domain_name
          target_domain_ary.push target_domain_hash[target_domain_name]
        else
          td = TargetDomain.new(target_domain_name, ie_cutoff, abs_len_cutoff, rel_len_cutoff)
          target_domain_hash[target_domain_name] = td
          target_domain_ary.push td
        end
      end
      target_domain_ary
    end

    def log(str = '')
      puts str
      @log << str << Constants::LINE_DELIMITER
    end

  end # class HmmscanMultiDomainExtractor

  class DomainData
    def initialize( name )
      if (name == nil) || name.size < 1
        error_msg = "domain name nil or empty"
        raise IOError, error_msg
      end
      @name = name
      @count = 0
      @i_e_value_min = 10000000.0
      @i_e_value_max = -1.0
      @i_e_value_sum = 0.0
      @len_min = 10000000
      @len_max = -1
      @len_sum = 0.0
    end

    def add( name, length, i_e_value)
      if name != @name
        error_msg = "domain names do not match"
        raise IOError, error_msg
      end

      if length < 0
        error_msg = "length cannot me negative"
        raise IOError, error_msg
      end
      if i_e_value < 0
        error_msg = "iE-value cannot me negative"
        raise IOError, error_msg
      end
      @count += 1
      @i_e_value_sum += i_e_value
      @len_sum += length
      if i_e_value > @i_e_value_max
        @i_e_value_max = i_e_value
      end
      if i_e_value < @i_e_value_min
        @i_e_value_min = i_e_value
      end
      if length > @len_max
        @len_max = length
      end
      if length < @len_min
        @len_min = length
      end
    end

    def avg_length
      if @count == 0
        return 0
      end
      @len_sum / @count
    end

    def avg_i_e_value
      if @count == 0
        return 0
      end
      @i_e_value_sum / @count
    end

    def to_str
      s = ''
      s << @name.rjust(16) + ': '
      s << @count.to_s.rjust(4) + '  '
      s << avg_length.round(1).to_s.rjust(6) + ' '
      s << @len_min.to_s.rjust(4) + ' -'
      s << @len_max.to_s.rjust(4) + '  '
      s << ("%.2E" % avg_i_e_value).rjust(9) + ' '
      s << ("%.2E" % @i_e_value_min).rjust(9) + ' -'
      s << ("%.2E" % @i_e_value_max).rjust(9) + ' '
      s
    end

    attr_reader :name, :count, :i_e_value_min, :i_e_value_max, :length_min, :length_max

  end

  class TargetDomain
    def initialize(name, i_e_value, abs_len, rel_len)
      if (name == nil) || name.size < 1
        error_msg = "target domain name nil or empty"
        raise IOError, error_msg
      end
      if rel_len > 1
        error_msg = name + ": target domain relative length is greater than 1"
        raise IOError, error_msg
      end
      if (abs_len <= 0) && (rel_len <= 0)
        error_msg = name + ": need to have either absolute length or relative length cutoff"
        raise IOError, error_msg
      end
      if (abs_len > 0) && (rel_len > 0)
        error_msg = name + ": cannot have both absolute length and relative length cutoff"
        raise IOError, error_msg
      end
      @name = name
      @i_e_value = i_e_value
      @abs_len = abs_len
      @rel_len = rel_len
    end

    def to_str
      s = @name.rjust(16) + ':'
      s << ' iE-cutoff: ' + ("%.2E" % @i_e_value).rjust(9)
      if @abs_len > 0
        s << ', abs len-cutoff: ' + @abs_len.to_s.rjust(4)
      end
      if @rel_len > 0
        s << ', rel len-cutoff: ' + @rel_len.to_s.rjust(4)
      end
      s
    end
    attr_reader :name, :i_e_value, :abs_len, :rel_len
  end

end # module Evoruby


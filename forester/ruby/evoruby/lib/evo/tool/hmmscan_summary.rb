#
# = lib/evo/tool/hmmscan_summary.rb - HmmscanSummary class
#
# Copyright::    Copyright (C) 2017 Christian M Zmasek
# License::      GNU Lesser General Public License (LGPL)

require 'set'
require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/io/parser/hmmscan_parser'

module Evoruby
  class HmmscanSummary

    PRG_NAME       = "hsp"
    PRG_VERSION    = "2.003"
    PRG_DESC       = "Summarize hmmscan output tables into simpler tables"
    PRG_DATE       = "170213"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    DELIMITER_OPTION              = "d"
    SPECIES_OPTION                = "s"
    I_E_VALUE_THRESHOLD_OPTION    = "ie"
    FS_E_VALUE_THRESHOLD_OPTION   = "pe"
    HMM_FOR_PROTEIN_OUTPUT        = "m"
    IGNORE_DUF_OPTION             = "i"
    PARSE_OUT_DESCRIPITION_OPTION = "a"
    HELP_OPTION_1                 = "help"
    HELP_OPTION_2                 = "h"

    USE_AVOID_HMMS = false
    AVOID_HHMS = [ "x", "y", "z" ]
    LIMIT_FOR_CLOSE_DOMAINS = 20 # Used for protein architecture summary

    def initialize
      @domain_counts = Hash.new
    end

    def run

      Util.print_program_information( PRG_NAME,
      PRG_VERSION,
      PRG_DESC,
      PRG_DATE,
      WWW,
      STDOUT )

      if ( ARGV == nil || ( ARGV.length < 1 )  )
        print_help
        exit( -1 )
      end

      begin
        cla = CommandLineArguments.new( ARGV )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      end

      if ( cla.is_option_set?( HELP_OPTION_1 ) ||
      cla.is_option_set?( HELP_OPTION_2 ) )
        print_help
        exit( 0 )
      end

      allowed_opts = Array.new
      allowed_opts.push( DELIMITER_OPTION )
      allowed_opts.push( I_E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( FS_E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( IGNORE_DUF_OPTION )
      allowed_opts.push( PARSE_OUT_DESCRIPITION_OPTION )
      allowed_opts.push( HMM_FOR_PROTEIN_OUTPUT )
      allowed_opts.push( SPECIES_OPTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
        "unknown option(s): " + disallowed,
        STDOUT )
      end

      inpath = cla.get_file_name( 0 )

      outpath = ""
      if ( cla.get_number_of_files == 1 )
        outpath = inpath + Constants::DOMAIN_TABLE_SUFFIX
      elsif ( cla.get_number_of_files == 2 )
        outpath = cla.get_file_name( 1 )
      else
        print_help
        exit( -1 )
      end

      column_delimiter = "\t"
      if ( cla.is_option_set?( DELIMITER_OPTION ) )
        begin
          column_delimiter = cla.get_option_value( DELIMITER_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      i_e_value_threshold = -1.0
      if ( cla.is_option_set?( I_E_VALUE_THRESHOLD_OPTION ) )
        begin
          i_e_value_threshold = cla.get_option_value_as_float( I_E_VALUE_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( i_e_value_threshold < 0.0 )
          Util.fatal_error( PRG_NAME, "attempt to use a negative i-E-value threshold", STDOUT )
        end
      end

      fs_e_value_threshold = -1.0
      if ( cla.is_option_set?( FS_E_VALUE_THRESHOLD_OPTION ) )
        begin
          fs_e_value_threshold = cla.get_option_value_as_float( FS_E_VALUE_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( fs_e_value_threshold < 0.0 )
          Util.fatal_error( PRG_NAME, "attempt to use a negative E-value threshold", STDOUT )
        end
      end

      hmm_for_protein_output = ""
      if ( cla.is_option_set?( HMM_FOR_PROTEIN_OUTPUT ) )
        begin
          hmm_for_protein_output = cla.get_option_value( HMM_FOR_PROTEIN_OUTPUT )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      species = "HUMAN"
      if ( cla.is_option_set?( SPECIES_OPTION ) )
        begin
          species = cla.get_option_value( SPECIES_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      ignore_dufs = false
      if ( cla.is_option_set?( IGNORE_DUF_OPTION ) )
        ignore_dufs = true
      end

      parse_descriptions = false
      if ( cla.is_option_set?( PARSE_OUT_DESCRIPITION_OPTION ) )
        parse_descriptions = true
      end

      puts()
      puts( "hmmscan outputfile  : " + inpath )
      puts( "outputfile          : " + outpath )

      if ( i_e_value_threshold >= 0.0 )
        puts( "i-E-value threshold : " + i_e_value_threshold.to_s )
      else
        puts( "i-E-value threshold : no threshold" )
      end
      if ( parse_descriptions )
        puts( "parse descriptions  : true" )
      else
        puts( "parse descriptions  : false" )
      end
      if ( ignore_dufs )
        puts( "ignore DUFs         : true" )
      else
        puts( "ignore DUFs         : false" )
      end
      if ( column_delimiter == "\t" )
        puts( "column delimiter    : TAB" )
      else
        puts( "column delimiter    : " + column_delimiter )
      end
      if !hmm_for_protein_output.empty?
        puts( "HMM for proteins    : " + hmm_for_protein_output )
        puts( "species             : " + species )
        if fs_e_value_threshold >= 0.0
          puts( "E-value threshold   : " + fs_e_value_threshold.to_s )
        else
          puts( "E-value threshold   : no threshold" )
        end
      end
      puts()

      begin
        parse( inpath,
        outpath,
        column_delimiter,
        i_e_value_threshold,
        ignore_dufs,
        parse_descriptions,
        fs_e_value_threshold,
        hmm_for_protein_output,
        species )
      rescue IOError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      end
      domain_counts = get_domain_counts()

      puts
      puts( "domain counts (considering potential i-E-value threshold and ignoring of DUFs):" )
      puts( "(number of different domains: " + domain_counts.length.to_s + ")" )
      puts
      puts( Util.draw_histogram( domain_counts, "#" ) )
      puts
      Util.print_message( PRG_NAME, "wrote: " + outpath )
      Util.print_message( PRG_NAME, "next step in standard analysis pipeline: d2f.rb")
      Util.print_message( PRG_NAME, 'OK' )
      puts

    end # def run

    private

    # raises ArgumentError, IOError
    def parse( inpath,
      outpath,
      column_delimiter,
      i_e_value_threshold,
      ignore_dufs,
      get_descriptions,
      fs_e_value_threshold,
      hmm_for_protein_output,
      species )

      Util.check_file_for_readability( inpath )
      Util.check_file_for_writability( outpath )

      hmmscan_parser = HmmscanParser.new( inpath )
      results = hmmscan_parser.parse

      outfile = File.open( outpath, "a" )

      query     = ""
      desc      = ""
      model     = ""
      env_from  = ""
      env_to    = ""
      i_e_value = ""

      hmmscan_results_per_protein = []

      prev_query = ""

      results.each do | r |
        model     = r.model
        desc      = r.desc
        query     = r.query
        i_e_value = r.i_e_value
        env_from  = r.env_from
        env_to    = r.env_to

        if ( ( i_e_value_threshold < 0.0 ) || ( i_e_value <= i_e_value_threshold ) ) &&
        ( !ignore_dufs || ( model !~ /^DUF\d+/ ) )
          count_model( model )
          outfile.print( query +
          column_delimiter )
          if ( get_descriptions )
            outfile.print( desc +
            column_delimiter )
          end
          outfile.print( model +
          column_delimiter +
          env_from.to_s +
          column_delimiter +
          env_to.to_s +
          column_delimiter +
          i_e_value.to_s )
          outfile.print( Constants::LINE_DELIMITER )
        end

        if !hmm_for_protein_output.empty?
          if  !prev_query.empty? && prev_query != query
            if !hmmscan_results_per_protein.empty?
              process_hmmscan_results_per_protein( hmmscan_results_per_protein,
              fs_e_value_threshold,
              hmm_for_protein_output,
              i_e_value_threshold,
              species )
            end
            hmmscan_results_per_protein.clear
          end
          prev_query = query

          if USE_AVOID_HMMS
            if !AVOID_HHMS.include? r.model
              hmmscan_results_per_protein << r
            end
          else
            hmmscan_results_per_protein << r
          end
        end
      end

      if !hmm_for_protein_output.empty? && !hmmscan_results_per_protein.empty?
        process_hmmscan_results_per_protein( hmmscan_results_per_protein,
        fs_e_value_threshold,
        hmm_for_protein_output,
        i_e_value_threshold,
        species )
      end

      outfile.flush()
      outfile.close()
    end # def parse

    def process_id( id )
      if  id =~ /(sp|tr)\|\S+\|(\S+)/
        id = $2
      end
      id
    end

    def count_model( model )
      if ( @domain_counts.has_key?( model ) )
        count = @domain_counts[ model ].to_i
        count += 1
        @domain_counts[ model ] = count
      else
        @domain_counts[ model ] = 1
      end
    end

    def process_hmmscan_results_per_protein( hmmscan_results_per_protein,
      fs_e_value_threshold,
      hmm_for_protein_output,
      i_e_value_threshold,
      species )

      dc = 0
      # filter according to i-Evalue threshold
      # abort if fs Evalue too high
      hmmscan_results_per_protein_filtered = []

      hmmscan_results_per_protein.each do | r |

        if r.model == hmm_for_protein_output
          if fs_e_value_threshold > 0.0 && r.fs_e_value > fs_e_value_threshold
            return
          end
        end
        if i_e_value_threshold <= 0 || r.i_e_value <= i_e_value_threshold
          hmmscan_results_per_protein_filtered << r
          if r.model == hmm_for_protein_output
            dc += 1
          end
        end
      end

      if dc == 0
        # passed on protein E-value, failed in per domain E-values
        return
      end

      hmmscan_results_per_protein_filtered.sort! { |r1,r2| r1.env_from <=> r2.env_from }

      own = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if r.model == hmm_for_protein_output
          own = r
        end
      end

      s = ""
      s << own.query + "\t"
      s << species + "\t"
      s << own.fs_e_value.to_s + "\t"
      s << own.qlen.to_s + "\t"
      s << dc.to_s + "\t"
      s << hmmscan_results_per_protein_filtered.length.to_s + "\t"
      hmmscan_results_per_protein_filtered.each do | r |
        s << r.model + " "
      end
      s << "\t"

      overview = make_overview( hmmscan_results_per_protein_filtered, hmm_for_protein_output )

      s << overview  + "\t"

      s << calc_linkers(  hmmscan_results_per_protein_filtered, hmm_for_protein_output )  + "\t"

      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if  prev_r != nil
          s << make_interdomain_sequence( r.env_from - prev_r.env_to - 1 )
        else
          s << make_interdomain_sequence( r.env_from, false )
        end
        s << r.model
        s << "["
        s << r.env_from.to_s << "-" << r.env_to.to_s
        s << "|ie=" << r.i_e_value.to_s
        s << "|ce=" << r.c_e_value.to_s
        s << "]"
        prev_r = r
      end
      s << make_interdomain_sequence( own.qlen - prev_r.env_to, false )
      puts s
    end

    def calc_linkers(  hmmscan_results_per_protein_filtered, hmm_for_protein_output )
      linkers = ""
      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if r.model == hmm_for_protein_output
          if  prev_r != nil
            linkers << ( r.env_from - prev_r.env_to - 1 ).to_s + " "
          end
          prev_r = r
        end
      end
      linkers
    end

    def get_domain_counts()
      return @domain_counts
    end

    def make_overview( hmmscan_results_per_protein_filtered, hmm_for_protein_output )
      overview = ""
      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if r.model == hmm_for_protein_output
          if prev_r == nil
            overview << hmm_for_protein_output
          else
            if  ( r.env_from - prev_r.env_to - 1 ) <= LIMIT_FOR_CLOSE_DOMAINS
              overview << "~" << hmm_for_protein_output
            else
              overview << "----" << hmm_for_protein_output
            end
          end
          prev_r = r
        end
      end
      overview
    end

    def make_interdomain_sequence( d, mark_short = true )
      s = ""
      d /= 20
      if d >= 10
        s << "----//----"
      elsif d >= 1
        d.times do
          s << "-"
        end
      elsif mark_short
        s << "~"
      end
      s
    end

    def print_help()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb [options] <hmmscan outputfile> [outputfile]" )
      puts()
      puts( "  options: -" + DELIMITER_OPTION + "=<s> : column delimiter for outputfile, default is TAB" )
      puts( "           -" + I_E_VALUE_THRESHOLD_OPTION  + "=<f>: i-E-value threshold, default is no threshold" )
      puts( "           -" + PARSE_OUT_DESCRIPITION_OPTION  + "     : parse query description (in addition to query name)" )
      puts( "           -" + IGNORE_DUF_OPTION  + "     : ignore DUFs" )
      puts( "           -" + HMM_FOR_PROTEIN_OUTPUT + "=<s> : HMM for protein architectures summary" )
      puts( "           -" + FS_E_VALUE_THRESHOLD_OPTION  + "=<f>: E-value threshold for full protein sequences, only for protein architectures summary" )
      puts( "           -" + SPECIES_OPTION + "=<s> : species for protein architectures summary" )
      puts()
      puts( "  [next step in standard analysis pipeline: d2f.rb]")
      puts()
      puts( "Examples:" )
      puts()
      puts( "  " + "hmmscan --max --domtblout P53_hmmscan_#{Constants::PFAM_V_FOR_EX}_10 -E 10 Pfam-A.hmm P53_ni.fasta" )
      puts()
      puts( "  " + PRG_NAME + ".rb P53_hmmscan_300_10" )
      puts()
    end

  end # class

end # module Evoruby
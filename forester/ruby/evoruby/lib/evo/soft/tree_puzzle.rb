#
# = lib/soft/tree_puzzle - TreePuzzle  class
#
# Copyright::  Copyright (C) 2009 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: tree_puzzle.rb,v 1.5 2009/10/08 22:44:54 cmzmasek Exp $
#
# last modified: 2009/10/06

require 'lib/evo/soft/resource_locations'
require 'lib/evo/util/util'

module Evoruby

    class TreePuzzle

        VERBOSE = false
        
        OUTDIST = 'outdist'
        OUTFILE = 'outfile'
        VERSION = '5.2'
        
        def initialize
            @tree_puzzle_home = Util.get_env_variable_value( ResourceLocations::TREEPUZZLE_HOME_ENV_VARIABLE )
            Util.check_file_for_readability( @tree_puzzle_home )
        end

        def run( alignment_file, model, rate_heterogeneity, number_of_seqs )
            Util.check_file_for_readability( alignment_file )

            input = alignment_file
            input << "\nk\nk"
            if number_of_seqs <= 257
                input << "\nk"
            end
            input << determine_model_option( model )
            input << determine_rate_heterogeneity_option( rate_heterogeneity )
            input << "\ny\n"

            if VERBOSE
                puts @tree_puzzle_home + " " + input
            end
            IO.popen( @tree_puzzle_home, 'r+' ) do |io|
                io.puts input
                io.close_write
                return io.read
            end
        end

        private

        # "Model of substitution" order for DQO TREE-PUZZLE 5.0:
        # Auto
        # m -> Dayhoff (Dayhoff et al. 1978)
        # m -> JTT (Jones et al. 1992)
        # m -> mtREV24 (Adachi-Hasegawa 1996)
        # m -> BLOSUM62 (Henikoff-Henikoff 92)
        # m -> VT (Mueller-Vingron 2000)
        # m -> WAG (Whelan-Goldman 2000)
        # m -> Auto
        def determine_model_option( model )
            cmd = nil
            if ( model == :pam )
                cmd = "\nm"
            elsif ( model == :jtt )
                cmd = "\nm\nm"
            elsif ( model == :mtrev24 )
                cmd = "\nm\nm\nm"
            elsif ( model == :blosum62 )
                cmd = "\nm\nm\nm\nm"
            elsif ( model == :vt )
                cmd = "\nm\nm\nm\nm\nm"
            elsif ( model == :wag )
                cmd = "\nm\nm\nm\nm\nm\nm"
            elsif ( model == :auto )
                cmd = ""
            else
                error_msg = "unknown model"
                raise ArgumentError, error_msg
            end
            cmd
        end


        # Model of rate heterogeneity:
        #    "8 Gamma distributed rates"
        #    "Two rates (1 invariable + 1 variable)"
        #    "Mixed (1 invariable + 8 Gamma rates)"
        #    otherwise: Uniform rate
        def determine_rate_heterogeneity_option( rates )
            opt = nil
            if ( rates == :gamma8 )
                opt = "\nw"
            elsif ( rates == :inv1_var1 )
                opt = "\nw\nw"
            elsif ( rates == :inv1_gamma8 )
                opt = "\nw\nw\nw"
            elsif ( rates == :uniform )
                opt = ""
            else
                error_msg = "unknown rate heterogeneity option"
                raise ArgumentError, error_msg
            end
            return opt
        end

    end # class TreePuzzle

end # module Evoruby

#!/usr/local/bin/ruby -w


require 'lib/evo/sequence/sequence'
require 'lib/evo/msa/msa'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/parser/fasta_parser'


module Evoruby

  input = ARGV[ 0 ]
  f = MsaFactory.new()

  IGNORE_SEQS_LACKING_GN = false

  msa = nil

  begin
    msa = f.create_msa_from_file( input, FastaParser.new() )
  rescue Exception => e
    puts "error: " + e.to_s
    exit
  end

  all_names = Set.new
  all_seqs_per_species = Hash.new
  gn_to_seqs = Hash.new
  unique_genes_msa = Msa.new
  longest_non_unique_genes_msa = Msa.new
  gn_re = /GN=(\S+)/
  fragment_re = /fragment/i
  species_re = /\[([A-Z0-9]{3,6})\]$/

  frag_counter = 0
  no_gn_counter = 0
  same_seq_counter = 0

  for i in 0 ... msa.get_number_of_seqs()
    seq = msa.get_sequence( i )
    name = seq.get_name
    if all_names.include?( name )
      puts "error: sequence name \"" + name + "\" is not unique (#" + i.to_s + ")"
      exit
    else
      all_names << name
    end

    if fragment_re.match( name )
      puts "ignored because fragment: " + name
      frag_counter += 1
      next
    end

    if species_re.match( name )
      s_match = species_re.match( name )
      species = s_match[1]

      unless all_seqs_per_species.include?( species )
        all_seqs_per_species[ species ] = Set.new
      end
      all_seqs = all_seqs_per_species[ species ]
      mol_seq = seq.get_sequence_as_string.upcase
      if all_seqs.include?( mol_seq )
        puts "ignored because identical sequence in same species: " + name
        same_seq_counter += 1
        next
      else
        all_seqs << mol_seq
      end
    else
      puts "no species for: " + name
    end

    gn_match = gn_re.match( name )
    if IGNORE_SEQS_LACKING_GN
      unless gn_match
        puts "ignored because no GN=: " + name
        no_gn_counter += 1
        next
      end
    else
      unless gn_match
        puts "no GN=: " + name
      end
    end

    gn =nil
    if gn_match
      gn = gn_match[1]
    else
      if IGNORE_SEQS_LACKING_GN
        puts "cannot be"
        exit
      end
      gn = name
    end

    unless gn_to_seqs.has_key?(gn)
      gn_to_seqs[gn] = Msa.new
    end
    gn_to_seqs[gn].add_sequence(seq)
  end

  puts "Sequences ignored because \"fragment\" in desc                : " + frag_counter.to_s
  if IGNORE_SEQS_LACKING_GN
    puts "Sequences ignored because no \"GN=\" in desc                  : " + no_gn_counter.to_s
  end
  puts "Sequences ignored because identical sequence in same species: " + same_seq_counter.to_s
  puts
  puts

  counter = 1
  gn_to_seqs.each_pair do |gene,seqs|
    if seqs.get_number_of_seqs > 1
      puts counter.to_s + ": " + gene
      puts seqs.to_fasta
      puts
      puts
      counter += 1
      longest = 0
      longest_seq = nil
      for j in 0 ... seqs.get_number_of_seqs()
        current = seqs.get_sequence( j )
        if current.get_length > longest
          longest =  current.get_length
          longest_seq = current
        end
      end
      longest_non_unique_genes_msa.add_sequence(longest_seq)
    else
      unique_genes_msa.add_sequence( seqs.get_sequence( 0 ) )
    end
  end
  w = FastaWriter.new
  w.write(unique_genes_msa, "seqs_from_unique_genes.fasta")
  w.write(longest_non_unique_genes_msa, "longest_seqs_from_nonunique_genes.fasta")
end

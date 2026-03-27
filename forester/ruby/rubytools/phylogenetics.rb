#!/usr/bin/env ruby

require 'set'
require 'optparse'

# Distance metrics for phylogenetic analysis
module DistanceMetrics
  PURINES = Set['A', 'G'].freeze
  PYRIMIDINES = Set['C', 'T', 'U'].freeze
  VALID_DNA = Set['A', 'C', 'G', 'T', 'U'].freeze
  VALID_AA = Set['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'].freeze
  AMBIGUITY_DNA = Set['N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'].freeze
  AMBIGUITY_AA = Set['X', 'B', 'Z', 'J'].freeze

  DNA_METRICS = %i[hamming jc k2p tn].freeze
  PROTEIN_METRICS = %i[p_distance poisson].freeze

  # Positions where both sequences have valid, unambiguous, non-gap characters
  def self.valid_dna_positions(seq1, seq2)
    valid_positions(seq1, seq2, VALID_DNA, AMBIGUITY_DNA)
  end

  def self.valid_aa_positions(seq1, seq2)
    valid_positions(seq1, seq2, VALID_AA, AMBIGUITY_AA)
  end

  # --- DNA Metrics ---

  # Simple proportional distance (gaps/ambiguities excluded)
  def self.hamming(seq1, seq2)
    pos = valid_dna_positions(seq1, seq2)
    return 0.0 if pos.empty?
    mismatches = pos.count { |i| seq1[i].upcase != seq2[i].upcase }
    mismatches.to_f / pos.length
  end

  # Jukes-Cantor: d = -3/4 * ln(1 - 4/3 * p)
  def self.jukes_cantor(seq1, seq2)
    p = hamming(seq1, seq2)
    return 0.0 if p == 0.0
    return Float::INFINITY if p >= 0.75
    -0.75 * Math.log(1.0 - (4.0 / 3.0) * p)
  end

  # Kimura 2-parameter: distinguishes transitions from transversions
  def self.kimura_2p(seq1, seq2)
    pos = valid_dna_positions(seq1, seq2)
    return 0.0 if pos.empty?
    ts, tv = count_ts_tv(seq1, seq2, pos)
    n = pos.length.to_f
    p = ts / n
    q = tv / n
    return 0.0 if p == 0.0 && q == 0.0
    a1 = 1.0 - 2.0 * p - q
    a2 = 1.0 - 2.0 * q
    return Float::INFINITY if a1 <= 0.0 || a2 <= 0.0
    -0.5 * Math.log(a1) - 0.25 * Math.log(a2)
  end

  # Tajima-Nei (1984): accounts for base frequency bias
  def self.tajima_nei(seq1, seq2)
    pos = valid_dna_positions(seq1, seq2)
    return 0.0 if pos.empty?

    freq = { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0 }
    pos.each do |i|
      freq[seq1[i].upcase] += 1
      freq[seq2[i].upcase] += 1
    end
    total = freq.values.sum.to_f
    freq.transform_values! { |v| v / total }

    p = pos.count { |i| seq1[i].upcase != seq2[i].upcase }.to_f / pos.length
    return 0.0 if p == 0.0

    b = 1.0 - freq.values.sum { |f| f * f }
    return Float::INFINITY if b == 0.0

    w = 1.0 - p / b
    return Float::INFINITY if w <= 0.0

    -b * Math.log(w)
  end

  # --- Protein Metrics ---

  # Proportion of differing amino acid positions
  def self.p_distance(seq1, seq2)
    pos = valid_aa_positions(seq1, seq2)
    return 0.0 if pos.empty?
    mismatches = pos.count { |i| seq1[i].upcase != seq2[i].upcase }
    mismatches.to_f / pos.length
  end

  # Poisson-corrected distance: d = -ln(1 - p)
  def self.poisson(seq1, seq2)
    p = p_distance(seq1, seq2)
    return 0.0 if p == 0.0
    return Float::INFINITY if p >= 0.95
    -Math.log(1.0 - p)
  end

  # --- Auto-detection ---

  def self.detect_type(sequences)
    sample = sequences.values.join.upcase.gsub(/[-.\sNX]/, '')
    return :dna if sample.empty?
    dna_count = sample.count('ACGTU')
    dna_count.to_f / sample.length > 0.9 ? :dna : :protein
  end

  # --- Internal helpers ---

  def self.valid_positions(seq1, seq2, valid_chars, ambiguity_chars)
    positions = []
    seq1.each_char.with_index do |c1, i|
      c2 = seq2[i]
      uc1 = c1.upcase
      uc2 = c2.upcase
      next if uc1 == '-' || uc2 == '-' || uc1 == '.' || uc2 == '.'
      next if ambiguity_chars.include?(uc1) || ambiguity_chars.include?(uc2)
      next unless valid_chars.include?(uc1) && valid_chars.include?(uc2)
      positions << i
    end
    positions
  end
  private_class_method :valid_positions

  def self.count_ts_tv(seq1, seq2, positions)
    ts = 0
    tv = 0
    positions.each do |i|
      c1 = seq1[i].upcase
      c2 = seq2[i].upcase
      next if c1 == c2
      if (PURINES.include?(c1) && PURINES.include?(c2)) ||
         (PYRIMIDINES.include?(c1) && PYRIMIDINES.include?(c2))
        ts += 1
      else
        tv += 1
      end
    end
    [ts, tv]
  end
  private_class_method :count_ts_tv
end

# Sequence file parser (FASTA, PHYLIP, stdin)
module SequenceParser
  def self.parse(source)
    content = if source == '-'
                $stdin.read
              else
                File.read(source)
              end

    lines = content.lines.map(&:chomp)
    return {} if lines.empty?

    # PHYLIP starts with a line like "5 9" or " 5 9"
    if lines.first.strip.match?(/^\d+\s+\d+\s*$/)
      parse_phylip(lines)
    else
      parse_fasta(lines)
    end
  end

  def self.parse_fasta(lines)
    sequences = {}
    current_name = nil
    current_seq = ""

    lines.each do |line|
      next if line.strip.empty?
      next if line.start_with?(';')

      if line.start_with?('>')
        sequences[current_name] = current_seq unless current_name.nil?
        current_name = line[1..].strip.split.first
        current_seq = ""
      else
        current_seq += line.strip
      end
    end
    sequences[current_name] = current_seq unless current_name.nil?
    sequences
  end

  def self.parse_phylip(lines)
    header = lines.first.strip.split
    ntaxa = header[0].to_i
    sequences = {}

    lines[1..].each do |line|
      next if line.strip.empty?
      break if sequences.length >= ntaxa

      parts = line.strip.split(/\s+/, 2)
      next unless parts.length == 2

      name = parts[0]
      seq = parts[1].gsub(/\s/, '')
      sequences[name] = seq
    end
    sequences
  end
end

# Node in a phylogenetic tree (supports multifurcations)
class TreeNode
  attr_accessor :name, :children, :distance, :bootstrap

  def initialize(name = nil, distance = 0.0)
    @name = name
    @children = []
    @distance = distance
    @bootstrap = nil
  end

  def leaf?
    @children.empty?
  end

  def taxa
    if leaf?
      [@name]
    else
      @children.flat_map(&:taxa)
    end
  end

  def to_newick
    if leaf?
      "#{escape_newick(@name)}:#{format_distance}"
    else
      inner = @children.map(&:to_newick).join(",")
      bs = @bootstrap ? @bootstrap.to_s : ""
      "(#{inner})#{bs}:#{format_distance}"
    end
  end

  def to_phyloxml(indent = "    ")
    xml = "#{indent}<clade>\n"
    xml += "#{indent}  <name>#{escape_xml(@name)}</name>\n" if @name
    xml += "#{indent}  <branch_length>#{format_distance}</branch_length>\n"
    if @bootstrap
      xml += "#{indent}  <confidence type=\"bootstrap\">#{@bootstrap}</confidence>\n"
    end
    @children.each { |c| xml += c.to_phyloxml(indent + "  ") }
    xml + "#{indent}</clade>\n"
  end

  def to_s(depth = 0)
    if leaf?
      "#{' ' * depth}#{@name} (d=#{format_distance})\n"
    else
      bs = @bootstrap ? " [#{@bootstrap}%]" : ""
      str = "#{' ' * depth}+#{bs}\n"
      @children.each { |c| str += c.to_s(depth + 2) }
      str
    end
  end

  private

  def format_distance
    d = @distance || 0.0
    d = 0.0 if d.respond_to?(:nan?) && d.nan?
    d = 999.0 if d.respond_to?(:infinite?) && d.infinite?
    d = [d, 0.0].max
    d.round(6).to_s
  end

  def escape_xml(str)
    str.to_s
       .gsub("&", "&amp;")
       .gsub("<", "&lt;")
       .gsub(">", "&gt;")
       .gsub("\"", "&quot;")
       .gsub("'", "&apos;")
  end

  def escape_newick(str)
    return "" if str.nil?
    if str.match?(/[\s,;:()\[\]']/)
      "'#{str.gsub("'", "''")}'"
    else
      str
    end
  end
end

# Colorful terminal progress bar
class ProgressBar
  BLOCK_CHARS = [" ", "\u258F", "\u258E", "\u258D", "\u258C", "\u258B", "\u258A", "\u2589", "\u2588"].freeze

  def initialize(total, label: "", width: 30, output: $stderr)
    @total = total
    @current = 0
    @label = label
    @width = width
    @output = output
    @start_time = Time.now
    @tty = output.respond_to?(:tty?) && output.tty?
    render
  end

  def advance(n = 1)
    @current += n
    render
  end

  def finish
    @current = @total
    render
    @output.puts
  end

  private

  def render
    return render_simple unless @tty

    pct = @current.to_f / @total
    filled_exact = pct * @width
    filled_full = filled_exact.floor
    frac_index = ((filled_exact - filled_full) * 8).round
    empty = @width - filled_full - (frac_index > 0 ? 1 : 0)

    # Color gradient: red (0%) → yellow (50%) → green (100%)
    r, g, b = gradient_color(pct)
    fg = "\e[38;2;#{r};#{g};#{b}m"
    bar_bg = "\e[48;2;40;40;40m"
    reset = "\e[0m"
    bold = "\e[1m"
    dim = "\e[2m"

    bar = "#{fg}#{bar_bg}"
    bar += BLOCK_CHARS[8] * filled_full
    bar += BLOCK_CHARS[frac_index] if frac_index > 0
    bar += " " * empty
    bar += reset

    elapsed = Time.now - @start_time
    eta_str = if @current > 0 && @current < @total
                remaining = elapsed / @current * (@total - @current)
                format_time(remaining)
              elsif @current >= @total
                format_time(elapsed)
              else
                "--:--"
              end

    status = @current >= @total ? "#{dim}elapsed#{reset}" : "#{dim}eta#{reset}"
    pct_color = "\e[38;2;#{r};#{g};#{b}m#{bold}"

    @output.print "\r  #{bold}#{@label}#{reset} #{bar} #{pct_color}#{(pct * 100).floor.to_s.rjust(3)}%#{reset} " \
                  "#{dim}#{@current}/#{@total}#{reset}  #{eta_str} #{status}  "
  end

  def render_simple
    if @current == 1
      @output.print "  #{@label} [#{@total}] "
    end
    @output.print "." if (@current % [(@total / 40.0).ceil, 1].max).zero?
  end

  def gradient_color(pct)
    if pct < 0.5
      t = pct * 2.0
      r = 255
      g = (100 + 155 * t).round
      b = 50
    else
      t = (pct - 0.5) * 2.0
      r = (255 - 200 * t).round
      g = 255
      b = (50 + 100 * t).round
    end
    [r.clamp(0, 255), g.clamp(0, 255), b.clamp(0, 255)]
  end

  def format_time(seconds)
    if seconds < 60
      "#{seconds.round(1)}s"
    elsif seconds < 3600
      "%d:%02d" % [seconds / 60, seconds % 60]
    else
      "%d:%02d:%02d" % [seconds / 3600, (seconds % 3600) / 60, seconds % 60]
    end
  end
end

# Phylogenetic tree inference engine
class PhylogeneticTree
  attr_reader :root, :distance_matrix, :taxa_names, :algorithm

  def initialize(sequences, algorithm: :upgma, metric: :hamming)
    @sequences = sequences
    @algorithm = algorithm.to_sym
    @metric = metric.to_sym
    @taxa_names = nil
    @distance_matrix = nil
    @rooted = (@algorithm == :upgma)
    @root = build_tree
  end

  def rooted?
    @rooted
  end

  def calculate_distance(seq1, seq2)
    case @metric
    when :hamming    then DistanceMetrics.hamming(seq1, seq2)
    when :jc         then DistanceMetrics.jukes_cantor(seq1, seq2)
    when :k2p        then DistanceMetrics.kimura_2p(seq1, seq2)
    when :tn         then DistanceMetrics.tajima_nei(seq1, seq2)
    when :p_distance then DistanceMetrics.p_distance(seq1, seq2)
    when :poisson    then DistanceMetrics.poisson(seq1, seq2)
    else raise "Unknown distance metric: #{@metric}"
    end
  end

  # --- Distance matrix ---

  def compute_distance_matrix
    @taxa_names = @sequences.keys
    n = @taxa_names.length
    @distance_matrix = Array.new(n) { Array.new(n, 0.0) }

    (0...n).each do |i|
      ((i + 1)...n).each do |j|
        d = calculate_distance(@sequences[@taxa_names[i]], @sequences[@taxa_names[j]])
        @distance_matrix[i][j] = d
        @distance_matrix[j][i] = d
      end
    end
  end

  def print_distance_matrix
    compute_distance_matrix if @distance_matrix.nil?
    max_len = @taxa_names.map(&:length).max
    col_w = [max_len, 10].max + 2

    print " " * (max_len + 2)
    @taxa_names.each { |n| printf "%-#{col_w}s", n }
    puts

    @taxa_names.each_with_index do |name, i|
      printf "%-#{max_len + 2}s", name
      @taxa_names.each_with_index do |_, j|
        d = @distance_matrix[i][j]
        if d == 0.0
          printf "%-#{col_w}s", "0.0"
        elsif d.infinite?
          printf "%-#{col_w}s", "Inf"
        else
          printf "%-#{col_w}.6f", d
        end
      end
      puts
    end
  end

  # --- Tree construction ---

  def build_tree
    compute_distance_matrix
    case @algorithm
    when :upgma then build_upgma
    when :nj    then build_nj
    else raise "Unknown algorithm: #{@algorithm}"
    end
  end

  def build_upgma
    matrix = @distance_matrix.map(&:dup)
    n = @taxa_names.length

    nodes = @taxa_names.map { |name| TreeNode.new(name, 0.0) }
    cluster_sizes = Array.new(n, 1)
    heights = Array.new(n, 0.0)
    active = (0...n).to_a

    while active.length > 1
      # Find minimum distance pair
      min_d = Float::INFINITY
      mi, mj = active[0], active[1]
      active.combination(2) do |i, j|
        if matrix[i][j] < min_d
          min_d = matrix[i][j]
          mi, mj = i, j
        end
      end

      # New node with correct branch lengths
      new_height = min_d / 2.0
      nodes[mi].distance = [new_height - heights[mi], 0.0].max
      nodes[mj].distance = [new_height - heights[mj], 0.0].max

      new_node = TreeNode.new
      new_node.children = [nodes[mi], nodes[mj]]

      # Update distances using UPGMA recurrence:
      # d(i∪j, k) = (|i|*d(i,k) + |j|*d(j,k)) / (|i|+|j|)
      remaining = active.reject { |x| x == mi || x == mj }
      ni = cluster_sizes[mi]
      nj = cluster_sizes[mj]
      remaining.each do |k|
        d = (ni * matrix[mi][k] + nj * matrix[mj][k]) / (ni + nj).to_f
        matrix[mi][k] = d
        matrix[k][mi] = d
      end

      # Reuse slot mi for merged cluster
      nodes[mi] = new_node
      cluster_sizes[mi] = ni + nj
      heights[mi] = new_height
      active = [mi] + remaining
    end

    nodes[active[0]]
  end

  def build_nj
    matrix = @distance_matrix.map(&:dup)
    n = @taxa_names.length

    # Replace infinite distances with a large finite value to avoid NaN
    max_finite = 0.0
    matrix.each { |row| row.each { |d| max_finite = d if d.finite? && d > max_finite } }
    cap = [max_finite * 10, 100.0].max
    matrix.each_with_index do |row, i|
      row.each_with_index do |d, j|
        matrix[i][j] = cap if !d.finite?
      end
    end

    nodes = @taxa_names.map { |name| TreeNode.new(name, 0.0) }
    active = (0...n).to_a

    while active.length > 2
      m = active.length

      # Net divergence for each node
      r = {}
      active.each do |i|
        r[i] = active.sum { |j| matrix[i][j] } / (m - 2).to_f
      end

      # Find minimum Q pair
      min_q = Float::INFINITY
      mi, mj = active[0], active[1]
      active.combination(2) do |i, j|
        q = matrix[i][j] - r[i] - r[j]
        if q < min_q
          min_q = q
          mi, mj = i, j
        end
      end

      # Branch lengths
      d_ij = matrix[mi][mj]
      branch_i = 0.5 * d_ij + 0.5 * (r[mi] - r[mj])
      branch_j = d_ij - branch_i
      nodes[mi].distance = branch_i.finite? ? [branch_i, 0.0].max : 0.0
      nodes[mj].distance = branch_j.finite? ? [branch_j, 0.0].max : 0.0

      new_node = TreeNode.new
      new_node.children = [nodes[mi], nodes[mj]]

      # Update distances: d(u,k) = (d(i,k) + d(j,k) - d(i,j)) / 2
      remaining = active.reject { |x| x == mi || x == mj }
      remaining.each do |k|
        d = (matrix[mi][k] + matrix[mj][k] - d_ij) / 2.0
        matrix[mi][k] = d
        matrix[k][mi] = d
      end

      nodes[mi] = new_node
      active = [mi] + remaining
    end

    # Final pair: split remaining distance equally
    i, j = active
    d = matrix[i][j] / 2.0
    d = 0.0 unless d.finite?
    nodes[i].distance = [d, 0.0].max
    nodes[j].distance = [d, 0.0].max

    final = TreeNode.new
    final.children = [nodes[i], nodes[j]]
    final
  end

  # --- Bootstrap ---

  def bootstrap(replicates)
    bar = ProgressBar.new(replicates, label: "Bootstrap")

    trees = []
    replicates.times do |i|
      resampled = resample_columns
      tree = PhylogeneticTree.new(resampled, algorithm: @algorithm, metric: @metric)
      trees << tree.root
      bar.advance
    end
    bar.finish

    annotate_bootstrap(root, trees, replicates)
    trees
  end

  def build_consensus(bootstrap_trees, total_replicates)
    all_taxa = root.taxa.sort
    total_set = all_taxa.to_set
    total_count = all_taxa.length

    # Count clade frequencies across all bootstrap trees
    clade_counts = Hash.new(0)
    bootstrap_trees.each do |bt|
      extract_clades(bt, total_count).each { |c| clade_counts[c] += 1 }
    end

    # Keep clades with >50% support
    supported = {}
    clade_counts.each do |clade, count|
      pct = ((count.to_f / total_replicates) * 100).round.to_i
      supported[clade] = pct if pct > 50
    end
    supported[total_set] = 100

    # Sort ascending by size so find() returns smallest match
    sorted = supported.keys.sort_by(&:length)

    # Create a node for each supported clade
    clade_nodes = {}
    sorted.each do |clade|
      node = TreeNode.new
      node.bootstrap = supported[clade]
      clade_nodes[clade] = node
    end

    # Assign each leaf to its smallest containing clade
    all_taxa.each do |taxon|
      parent_clade = sorted.find { |c| c.include?(taxon) }
      clade_nodes[parent_clade].children << TreeNode.new(taxon) if parent_clade
    end

    # Assign each clade to its smallest proper containing clade
    sorted.each do |clade|
      parent_clade = sorted.find { |c| c.proper_superset?(clade) }
      clade_nodes[parent_clade].children << clade_nodes[clade] if parent_clade
    end

    clade_nodes[total_set]
  end

  # --- Rooting ---

  # Re-root tree using outgroup taxon
  def reroot_by_outgroup!(outgroup_name)
    outgroup_leaf = find_leaf(@root, outgroup_name)
    raise "Outgroup '#{outgroup_name}' not found" unless outgroup_leaf

    adj = build_adjacency_list(@root)

    # Outgroup's only neighbor is its parent
    parent_node, edge_dist = adj[outgroup_leaf].first

    new_root = TreeNode.new

    # Outgroup side
    og = TreeNode.new(outgroup_name, edge_dist)
    new_root.children << og

    # Ingroup side: rebuild from parent excluding outgroup direction
    ingroup = build_rooted_from_adj(adj, parent_node, outgroup_leaf)
    ingroup.distance = 0.0
    new_root.children << ingroup

    @root = collapse_degree_one(new_root)
    @rooted = true
  end

  # Re-root tree at midpoint of longest path
  def reroot_by_midpoint!
    adj = build_adjacency_list(@root)
    leaves = collect_leaf_nodes(@root)

    # Find diameter: two most distant leaves
    max_dist = 0.0
    far_a = far_b = leaves.first

    # Two-pass optimization: BFS from arbitrary leaf to find farthest,
    # then BFS from that to find the true farthest pair
    dists_from_first = bfs_distances(adj, leaves.first)
    far_a = leaves.max_by { |l| dists_from_first[l] || 0.0 }

    dists_from_a = bfs_distances(adj, far_a)
    far_b = leaves.max_by { |l| dists_from_a[l] || 0.0 }
    max_dist = dists_from_a[far_b] || 0.0

    if max_dist == 0.0
      @rooted = true
      return
    end

    # Find path between the two most distant leaves
    path = find_path_in_adj(adj, far_a, far_b)
    return unless path

    # Walk along path to midpoint (diameter / 2 from far_a)
    target = max_dist / 2.0
    accumulated = 0.0

    (0...(path.length - 1)).each do |i|
      from_node = path[i]
      to_node = path[i + 1]
      edge_dist = adj[from_node].find { |n, _| n == to_node }[1]

      if accumulated + edge_dist >= target - 1e-10
        # Split this edge at the midpoint
        dist_from_a_side = target - accumulated
        dist_from_b_side = edge_dist - dist_from_a_side

        # Clamp to avoid negative from floating point
        dist_from_a_side = [dist_from_a_side, 0.0].max
        dist_from_b_side = [dist_from_b_side, 0.0].max

        new_root = TreeNode.new

        # Side toward far_a
        side_a = build_rooted_from_adj(adj, from_node, to_node)
        side_a.distance = dist_from_a_side
        new_root.children << side_a

        # Side toward far_b
        side_b = build_rooted_from_adj(adj, to_node, from_node)
        side_b.distance = dist_from_b_side
        new_root.children << side_b

        @root = collapse_degree_one(new_root)
        @rooted = true
        return
      end

      accumulated += edge_dist
    end

    @rooted = true
  end

  # --- Output ---

  def to_phyloxml
    rooted_str = rooted? ? "true" : "false"
    xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    xml += "<phyloxml xmlns=\"http://www.phyloxml.org\" "
    xml += "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
    xml += "xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\">\n"
    xml += "  <phylogeny rooted=\"#{rooted_str}\">\n"
    xml += root.to_phyloxml("    ")
    xml += "  </phylogeny>\n"
    xml += "</phyloxml>\n"
    xml
  end

  private

  def resample_columns
    seq_length = @sequences.values.first.length
    columns = Array.new(seq_length) { rand(seq_length) }
    @sequences.transform_values { |seq| columns.map { |c| seq[c] }.join }
  end

  def annotate_bootstrap(node, bootstrap_trees, total)
    return if node.leaf?

    taxa_set = node.taxa.to_set
    count = bootstrap_trees.count { |bt| has_clade?(bt, taxa_set) }
    node.bootstrap = ((count.to_f / total) * 100).round.to_i

    node.children.each { |c| annotate_bootstrap(c, bootstrap_trees, total) }
  end

  def has_clade?(node, taxa_set)
    return false if node.leaf?
    return true if node.taxa.to_set == taxa_set
    node.children.any? { |c| has_clade?(c, taxa_set) }
  end

  def extract_clades(node, total_count)
    return [] if node.leaf?
    clades = []
    taxa = node.taxa.to_set
    clades << taxa if taxa.length > 1 && taxa.length < total_count
    node.children.each { |c| clades += extract_clades(c, total_count) }
    clades
  end

  # --- Rooting helpers ---

  def find_leaf(node, name)
    return node if node.leaf? && node.name == name
    node.children.each do |child|
      result = find_leaf(child, name)
      return result if result
    end
    nil
  end

  def collect_leaf_nodes(node)
    return [node] if node.leaf?
    node.children.flat_map { |c| collect_leaf_nodes(c) }
  end

  # Convert tree to undirected adjacency list (edge weights = branch lengths)
  def build_adjacency_list(node)
    adj = {}
    _adj_recurse(node, adj)
    adj
  end

  def _adj_recurse(node, adj)
    adj[node] ||= []
    node.children.each do |child|
      adj[node] << [child, child.distance]
      adj[child] ||= []
      adj[child] << [node, child.distance]
      _adj_recurse(child, adj)
    end
  end

  # Build a rooted tree from adjacency list starting at a node
  def build_rooted_from_adj(adj, current, came_from)
    node = TreeNode.new(current.leaf? ? current.name : nil)
    node.bootstrap = current.bootstrap
    node.distance = 0.0

    adj[current].each do |neighbor, dist|
      next if neighbor == came_from
      child = build_rooted_from_adj(adj, neighbor, current)
      child.distance = dist
      node.children << child
    end

    node
  end

  # BFS to compute distances from a start node to all others
  def bfs_distances(adj, start)
    distances = { start => 0.0 }
    queue = [start]
    until queue.empty?
      current = queue.shift
      adj[current].each do |neighbor, dist|
        next if distances.key?(neighbor)
        distances[neighbor] = distances[current] + dist
        queue << neighbor
      end
    end
    distances
  end

  # Find path between two nodes in adjacency list (DFS)
  def find_path_in_adj(adj, from, to, visited = Set.new)
    return [from] if from == to
    visited.add(from)
    adj[from].each do |neighbor, _|
      next if visited.include?(neighbor)
      sub_path = find_path_in_adj(adj, neighbor, to, visited)
      return [from] + sub_path if sub_path
    end
    nil
  end

  # Collapse internal nodes with exactly one child
  def collapse_degree_one(node)
    node.children.map! { |c| collapse_degree_one(c) }

    if !node.leaf? && node.children.length == 1
      child = node.children[0]
      child.distance += node.distance
      child.bootstrap ||= node.bootstrap
      return child
    end

    node
  end
end

# Command-line interface
class PhyloCLI
  METRIC_ALIASES = {
    'hamming' => :hamming,
    'jc' => :jc, 'jukes-cantor' => :jc, 'jukes_cantor' => :jc,
    'k2p' => :k2p, 'kimura' => :k2p, 'kimura-2p' => :k2p, 'kimura_2p' => :k2p,
    'tn' => :tn, 'tajima-nei' => :tn, 'tajima_nei' => :tn,
    'p_distance' => :p_distance, 'p-distance' => :p_distance, 'pdist' => :p_distance,
    'poisson' => :poisson
  }.freeze

  def initialize
    @options = {
      algorithm: :upgma,
      metric: nil,
      format: :both,
      output: nil,
      demo: false,
      bootstrap: nil,
      consensus: false,
      matrix: false,
      seed: nil,
      type: nil,
      root: nil
    }
  end

  def parse_args
    OptionParser.new do |opts|
      opts.banner = "Usage: phylogenetics.rb [options]"
      opts.separator ""
      opts.separator "Input:"

      opts.on("-f", "--file FILE", "Input file (FASTA or PHYLIP; '-' for stdin)") do |f|
        @options[:file] = f
      end

      opts.on("-d", "--demo", "Run with built-in demo sequences") do
        @options[:demo] = true
      end

      opts.on("-t", "--type TYPE", "Sequence type: dna or protein (default: auto)") do |t|
        type = t.downcase.to_sym
        abort "Error: Type must be 'dna' or 'protein'" unless %i[dna protein].include?(type)
        @options[:type] = type
      end

      opts.separator ""
      opts.separator "Analysis:"

      opts.on("-a", "--algorithm ALG", "Tree algorithm: upgma or nj (default: upgma)") do |a|
        algo = a.downcase.to_sym
        abort "Error: Algorithm must be 'upgma' or 'nj'" unless %i[upgma nj].include?(algo)
        @options[:algorithm] = algo
      end

      opts.on("-m", "--metric METRIC",
              "Distance metric (default: jc for DNA, poisson for protein)",
              "  DNA:     hamming, jc, k2p, tn",
              "  Protein: p_distance, poisson") do |m|
        metric = METRIC_ALIASES[m.downcase]
        abort "Error: Unknown metric '#{m}'" unless metric
        @options[:metric] = metric
      end

      opts.on("-b", "--bootstrap NUM", Integer, "Number of bootstrap replicates") do |b|
        abort "Error: Bootstrap replicates must be >= 1" if b < 1
        @options[:bootstrap] = b
      end

      opts.on("-c", "--consensus", "Build majority-rule consensus tree (requires -b)") do
        @options[:consensus] = true
      end

      opts.on("-r", "--root METHOD",
              "Root NJ tree: midpoint or outgroup:NAME") do |r|
        if r.downcase == 'midpoint'
          @options[:root] = :midpoint
        elsif r.downcase.start_with?('outgroup:')
          name = r.split(':', 2)[1]
          abort "Error: Outgroup name required (e.g. --root outgroup:Fish)" if name.nil? || name.empty?
          @options[:root] = [:outgroup, name]
        else
          abort "Error: Root method must be 'midpoint' or 'outgroup:NAME'"
        end
      end

      opts.on("-s", "--seed NUM", Integer, "Random seed for reproducible bootstrap") do |s|
        @options[:seed] = s
      end

      opts.separator ""
      opts.separator "Output:"

      opts.on("-F", "--format FORMAT", "Output format: newick, tree, phyloxml, or both (default: both)") do |f|
        fmt = f.downcase.to_sym
        abort "Error: Format must be 'newick', 'tree', 'phyloxml', or 'both'" unless %i[newick tree phyloxml both].include?(fmt)
        @options[:format] = fmt
      end

      opts.on("-o", "--output FILE", "Save tree to file") do |o|
        @options[:output] = o
      end

      opts.on("-M", "--matrix", "Print the pairwise distance matrix") do
        @options[:matrix] = true
      end

      opts.on("-h", "--help", "Show this help message") do
        puts opts
        exit 0
      end
    end.parse!
  end

  def run
    parse_args

    abort "Error: --consensus requires --bootstrap" if @options[:consensus] && !@options[:bootstrap]
    if @options[:root] && @options[:algorithm] != :nj
      abort "Error: --root is only applicable to NJ trees (use -a nj)"
    end
    srand(@options[:seed]) if @options[:seed]

    sequences = load_sequences
    validate_sequences(sequences)

    seq_type = @options[:type] || DistanceMetrics.detect_type(sequences)
    metric = resolve_metric(seq_type)
    seq_length = sequences.values.first.length

    puts "#{sequences.length} taxa, #{seq_length} sites (#{seq_type.upcase})"
    puts "Algorithm: #{@options[:algorithm].to_s.upcase} | Metric: #{metric.to_s.upcase}"
    puts "=" * 60

    tree = PhylogeneticTree.new(sequences, algorithm: @options[:algorithm], metric: metric)

    # Root NJ tree if requested
    if @options[:root]
      case @options[:root]
      when :midpoint
        tree.reroot_by_midpoint!
        $stderr.puts "Rooted by midpoint"
      when Array
        outgroup_name = @options[:root][1]
        available = sequences.keys.join(', ')
        begin
          tree.reroot_by_outgroup!(outgroup_name)
        rescue RuntimeError => e
          abort "Error: #{e.message}. Available taxa: #{available}"
        end
        $stderr.puts "Rooted by outgroup: #{outgroup_name}"
      end
    end

    if @options[:matrix]
      puts "\nDistance Matrix:"
      tree.print_distance_matrix
    end

    bootstrap_trees = nil
    if @options[:bootstrap]
      bootstrap_trees = tree.bootstrap(@options[:bootstrap])
    end

    puts
    output_tree(tree)

    if @options[:consensus] && bootstrap_trees
      consensus = tree.build_consensus(bootstrap_trees, @options[:bootstrap])
      puts "\n#{"=" * 60}"
      puts "Majority-Rule Consensus Tree"
      puts "=" * 60
      puts "Newick:"
      puts consensus.to_newick
      puts "\nTree:"
      puts consensus.to_s
    end

    save_output(tree)
  end

  private

  def load_sequences
    if @options[:demo]
      {
        "Human"      => "ATGGCGCCC",
        "Chimpanzee" => "ATGGCGCCC",
        "Mouse"      => "ATGACGCCC",
        "Chicken"    => "ATGACGTAT",
        "Fish"       => "ATGAGATAT"
      }
    elsif @options[:file]
      unless @options[:file] == '-' || File.exist?(@options[:file])
        abort "Error: File not found: #{@options[:file]}"
      end
      SequenceParser.parse(@options[:file])
    else
      abort "Error: Provide a file with -f, use -d for demo, or pipe via -f -\nRun with -h for help"
    end
  end

  def validate_sequences(sequences)
    abort "Error: No sequences found" if sequences.empty?
    abort "Error: Need at least 2 sequences" if sequences.length < 2

    lengths = sequences.values.map(&:length).uniq
    return if lengths.length == 1

    abort "Error: Sequences must be aligned (same length). Found lengths: #{lengths.sort.join(', ')}"
  end

  def resolve_metric(seq_type)
    if @options[:metric]
      m = @options[:metric]
      if seq_type == :dna && DistanceMetrics::PROTEIN_METRICS.include?(m)
        $stderr.puts "Warning: Using protein metric '#{m}' on DNA sequences"
      elsif seq_type == :protein && DistanceMetrics::DNA_METRICS.include?(m)
        $stderr.puts "Warning: Using DNA metric '#{m}' on protein sequences"
      end
      m
    else
      seq_type == :protein ? :poisson : :jc
    end
  end

  def output_tree(tree)
    case @options[:format]
    when :newick
      puts tree.root.to_newick
    when :tree
      puts tree.root.to_s
    when :phyloxml
      puts tree.to_phyloxml
    when :both
      puts "Newick:"
      puts tree.root.to_newick
      puts "\nTree:"
      puts tree.root.to_s
    end
  end

  def save_output(tree)
    return unless @options[:output]

    content = case @options[:format]
              when :phyloxml then tree.to_phyloxml
              else tree.root.to_newick
              end
    File.write(@options[:output], content)
    $stderr.puts "Tree saved to: #{@options[:output]}"
  end
end

if __FILE__ == $0
  PhyloCLI.new.run
end
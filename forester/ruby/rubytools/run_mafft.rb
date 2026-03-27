#!/usr/bin/env ruby
# frozen_string_literal: true

VERSION = '1.0.0'

require 'optparse'
require 'fileutils'
require 'shellwords'

# ── Progress bar ──────────────────────────────────────────────────────────────

class ProgressBar
  BAR_WIDTH = 40

  def initialize(total)
    @total   = total
    @current = 0
    @start   = Time.now
    print "\n"
  end

  def increment(label = '')
    @current += 1
    render(label)
  end

  def finish
    render('done', final: true)
    print "\n\n"
  end

  def elapsed
    Time.now - @start
  end

  private

  def render(label, final: false)
    pct     = @total.zero? ? 1.0 : @current.to_f / @total
    filled  = (pct * BAR_WIDTH).round
    bar     = '#' * filled + '-' * (BAR_WIDTH - filled)
    elapsed = Time.now - @start
    eta     = (pct > 0 && !final) ? (elapsed / pct * (1 - pct)).round : 0
    eta_str = final ? '      ' : "ETA #{eta}s"
    label   = label.length > 30 ? "#{label[0..27]}..." : label.ljust(30)

    print "\r  [#{bar}] #{@current.to_s.rjust(@total.to_s.length)}/#{@total}  " \
          "#{(pct * 100).round(1).to_s.rjust(5)}%  #{eta_str}  #{label}"
    $stdout.flush
  end
end

# ── Helpers ───────────────────────────────────────────────────────────────────

def format_elapsed(seconds)
  seconds = seconds.round
  return "#{seconds}s" if seconds < 60

  mins, secs = seconds.divmod(60)
  return "#{mins}m #{secs}s" if mins < 60

  hrs, mins = mins.divmod(60)
  "#{hrs}h #{mins}m #{secs}s"
end

# Parse a FASTA file into an array of {id:, seq:} hashes.
def parse_fasta(path)
  records     = []
  current_id  = nil
  current_seq = +''
  File.foreach(path) do |line|
    line.chomp!
    if line.start_with?('>')
      records << { id: current_id, seq: current_seq } if current_id
      current_id  = line[1..].split.first
      current_seq = +''
    elsif current_id
      current_seq << line
    end
  end
  records << { id: current_id, seq: current_seq } if current_id
  records
end

# Validate a single FASTA file.  Returns an array of human-readable issue strings.
# label is a short descriptor used in messages, e.g. "aligned", "unaligned", "input".
def validate_fasta(path, label)
  issues  = []
  records = parse_fasta(path)

  if records.empty?
    issues << "#{label} file contains no sequences"
    return issues
  end

  # Duplicate sequence IDs
  id_counts = records.map { |r| r[:id] }.each_with_object(Hash.new(0)) { |id, h| h[id] += 1 }
  dups      = id_counts.select { |_, n| n > 1 }.keys
  issues << "#{label} file has duplicate ID(s): #{dups.join(', ')}" if dups.any?

  # Empty or all-gap/space sequences
  records.each do |r|
    residues = r[:seq].gsub(/[-.\s]/, '')
    if residues.empty?
      desc = r[:seq].empty? ? 'has no residues' : 'contains only gap/space characters'
      issues << "#{label} file: sequence '#{r[:id]}' #{desc}"
    end
  end

  issues
end

# Add-mode check: verify all sequences in the aligned file have equal length,
# i.e. it is actually an alignment and not a bag of unaligned sequences.
def check_equal_lengths(path, label)
  records = parse_fasta(path)
  return [] if records.size < 2   # single sequence is trivially "aligned"

  lengths = records.map { |r| r[:seq].length }
  return [] if lengths.uniq.size == 1

  min, max = lengths.minmax
  ["#{label} file is not a valid alignment: sequence lengths range from #{min} to #{max}"]
end

# Convert an array of {id:, seq:} records to PHYLIP sequential format.
# PHYLIP requires a fixed 10-character name field; longer names are truncated.
# Returns an array of warning strings (empty when everything is fine).
def write_phylip(records, path)
  warnings = []
  ntax     = records.size
  nchar    = records.first[:seq].length

  phylip_names = records.map do |r|
    warnings << "  PHYLIP: ID '#{r[:id]}' truncated to '#{r[:id][0, 10]}'" if r[:id].length > 10
    r[:id][0, 10].ljust(10)
  end

  dup_names = phylip_names.map(&:strip)
                          .each_with_object(Hash.new(0)) { |n, h| h[n] += 1 }
                          .select { |_, c| c > 1 }.keys
  warnings << "  PHYLIP: truncated names are not unique: #{dup_names.join(', ')}" if dup_names.any?

  File.open(path, 'w') do |f|
    f.puts " #{ntax} #{nchar}"
    records.each_with_index { |r, i| f.puts "#{phylip_names[i]}#{r[:seq]}" }
  end

  warnings
end

# Compute basic MSA quality metrics from an array of {id:, seq:} records.
# Returns a hash of metrics, or nil if the array is empty.
def compute_msa_metrics(records)
  return nil if records.empty?

  ntax  = records.size
  nchar = records.first[:seq].length

  gap_chars = /[-\.]/

  # Per-sequence gap percentages
  seq_gap_pcts = records.map do |r|
    { id: r[:id], gap_pct: (r[:seq].count('-') + r[:seq].count('.')) * 100.0 / nchar }
  end

  total_gaps      = records.sum { |r| r[:seq].count('-') + r[:seq].count('.') }
  overall_gap_pct = total_gaps * 100.0 / (ntax * nchar)
  mean_gap_pct    = seq_gap_pcts.sum { |s| s[:gap_pct] } / ntax
  most_gapped     = seq_gap_pcts.max_by { |s| s[:gap_pct] }
  least_gapped    = seq_gap_pcts.min_by { |s| s[:gap_pct] }

  # Column-wise analysis via transposition
  columns        = records.map { |r| r[:seq].chars }.transpose
  all_gap_cols   = columns.count { |col| col.all? { |c| c =~ gap_chars } }
  conserved_cols = columns.count do |col|
    non_gap = col.reject { |c| c =~ gap_chars }
    non_gap.size == ntax && non_gap.uniq.size == 1
  end

  # Columns where gap fraction exceeds a threshold
  cols_gap50  = columns.count { |col| col.count { |c| c =~ gap_chars } * 100.0 / ntax > 50 }
  cols_gap90  = columns.count { |col| col.count { |c| c =~ gap_chars } * 100.0 / ntax > 90 }

  {
    ntax:              ntax,
    nchar:             nchar,
    gap_pct:           overall_gap_pct.round(2),
    mean_seq_gap_pct:  mean_gap_pct.round(2),
    max_gap_pct:       most_gapped[:gap_pct].round(2),
    max_gap_seq:       most_gapped[:id],
    min_gap_pct:       least_gapped[:gap_pct].round(2),
    min_gap_seq:       least_gapped[:id],
    conserved_cols:    conserved_cols,
    all_gap_cols:      all_gap_cols,
    cols_gt50pct_gaps: cols_gap50,
    cols_gt50pct_gaps_pct: (cols_gap50 * 100.0 / nchar).round(2),
    cols_gt90pct_gaps: cols_gap90,
    cols_gt90pct_gaps_pct: (cols_gap90 * 100.0 / nchar).round(2)
  }
end

# Add-mode check: warn if any ID from the unaligned file already exists in the
# aligned file — MAFFT --add would silently produce a corrupted alignment.
def check_id_overlap(aligned_path, unaligned_path)
  aligned_ids   = parse_fasta(aligned_path).map  { |r| r[:id] }
  unaligned_ids = parse_fasta(unaligned_path).map { |r| r[:id] }
  overlap = aligned_ids & unaligned_ids
  return [] if overlap.empty?

  ["ID(s) already present in the aligned file: #{overlap.join(', ')}"]
end

def prompt(msg, default = nil)
  if default
    print "  #{msg} [#{default}]: "
  else
    print "  #{msg}: "
  end
  $stdout.flush
  input = $stdin.gets.chomp
  input.empty? ? default : input
end

# ── Option parsing ────────────────────────────────────────────────────────────

options = {
  input_dir:      nil,
  input_suffix:   '.fasta',
  output_dir:     nil,
  output_suffix:  '_mafft.fasta',
  mafft_opts:     '--auto',
  threads:        nil,
  overwrite:      false,
  dry_run:        false,
  add_mode:       false,
  phylip_suffix:  nil,
  msa_quality:    nil,
  aligned_suffix: '_aln.fasta',
  log_file_path:  nil
}

interactive    = ARGV.empty?       # capture before parse! mutates ARGV
original_argv  = ARGV.dup          # preserve for log reproduction

parser = OptionParser.new do |o|
  o.banner = <<~BANNER

    MAFFT Batch Aligner v#{VERSION}

    Usage: ruby run_mafft.rb [options]

    Standard mode:
      Runs MAFFT on every FASTA file matching --input-suffix and writes
      aligned output files to --output-dir.

    Add mode (--add-mode):
      Uses MAFFT's --add --keeplength to add unaligned sequences into an
      existing alignment, preserving its length. Files are matched by stem:
        <stem><aligned-suffix>   already-aligned reference
        <stem><input-suffix>     unaligned sequences to add
        <stem><output-suffix>    output

  BANNER

  o.separator '  Common options:'
  o.on('-i', '--input-dir DIR',        'Input directory (required)')               { |v| options[:input_dir]      = v }
  o.on('-o', '--output-dir DIR',       'Output directory (required)')              { |v| options[:output_dir]     = v }
  o.on('-x', '--output-suffix EXT',    'Output file suffix     [_mafft.fasta]')    { |v| options[:output_suffix]  = v }
  o.on('-m', '--mafft-opts "OPTS"',    'MAFFT options string   [--auto]')          { |v| options[:mafft_opts]     = v }
  o.on('-t', '--threads N',    Integer, 'Threads passed to MAFFT')                 { |v| options[:threads]        = v }
  o.on('-f', '--[no-]overwrite',       'Overwrite existing output files [false]')  { |v| options[:overwrite]      = v }
  o.on('-n', '--dry-run',              'Print commands without executing them')     {     options[:dry_run]          = true }
  o.on('-p', '--phylip-suffix EXT',    'Also write PHYLIP sequential output')       { |v| options[:phylip_suffix]  = v }
  o.on('-q', '--msa-quality FILE',     'Write MSA quality metrics TSV to FILE')     { |v| options[:msa_quality]    = v }

  o.separator ''
  o.separator '  Standard mode:'
  o.on('-s', '--input-suffix EXT',     'Input file suffix      [.fasta]')          { |v| options[:input_suffix]   = v }

  o.separator ''
  o.separator '  Add mode:'
  o.on('-A', '--add-mode',             'Enable --add --keeplength mode')           {     options[:add_mode]        = true }
  o.on('-a', '--aligned-suffix EXT',   'Already-aligned suffix [_aln.fasta]')      { |v| options[:aligned_suffix] = v }

  o.separator ''
  o.on('-l', '--log-file FILE',        'Write detailed run log to FILE')           { |v| options[:log_file_path] = v }
  o.on('-v', '--version',              'Print version and exit')                   { puts "MAFFT Batch Aligner v#{VERSION}"; exit }
  o.on('-h', '--help',                 'Show this help')                           { puts o; exit }
end

parser.parse!

# ── Banner ────────────────────────────────────────────────────────────────────

puts
puts '╔══════════════════════════════════════════╗'
puts "║       MAFFT Batch Aligner v#{VERSION.ljust(13)}║"
puts '╚══════════════════════════════════════════╝'
puts
$stdout.flush

# ── Interactive prompts (TTY only, no CLI args given) ─────────────────────────

if $stdin.tty? && interactive
  options[:input_dir]     = prompt('Input directory')
  options[:output_dir]    = prompt('Output directory')
  options[:output_suffix] = prompt('Output suffix',    options[:output_suffix])
  options[:mafft_opts]    = prompt('MAFFT options',    options[:mafft_opts])

  add = prompt('Enable add mode? (y/n)', 'n')
  if add.downcase.start_with?('y')
    options[:add_mode]       = true
    options[:aligned_suffix] = prompt('Aligned suffix',   options[:aligned_suffix])
    options[:input_suffix]   = prompt('Unaligned suffix', options[:input_suffix])
  else
    options[:input_suffix]   = prompt('Input suffix',     options[:input_suffix])
  end
end

# ── Validate required options ─────────────────────────────────────────────────

missing = []
missing << '--input-dir'  unless options[:input_dir]
missing << '--output-dir' unless options[:output_dir]
unless missing.empty?
  abort "\n  ERROR: required options missing: #{missing.join(', ')}\n" \
        "  Run with --help for usage.\n"
end

abort "\n  ERROR: input directory '#{options[:input_dir]}' not found.\n" \
  unless Dir.exist?(options[:input_dir])

FileUtils.mkdir_p(options[:output_dir]) unless options[:dry_run]

# ── Suffix collision checks (before any file discovery) ───────────────────────

if options[:add_mode]
  if options[:aligned_suffix] == options[:output_suffix]
    abort "\n  ERROR: aligned suffix and output suffix are identical ('#{options[:aligned_suffix]}').\n" \
          "  This would overwrite your aligned reference files. Use -x to set a different output suffix.\n"
  end
  if options[:input_suffix] == options[:output_suffix]
    abort "\n  ERROR: unaligned suffix and output suffix are identical ('#{options[:input_suffix]}').\n" \
          "  Use -x to set a different output suffix.\n"
  end
else
  if options[:input_suffix] == options[:output_suffix]
    abort "\n  ERROR: input suffix and output suffix are identical ('#{options[:input_suffix]}').\n" \
          "  This would overwrite your source files. Use -x to set a different output suffix.\n"
  end
end

# Phylip suffix must not clash with FASTA output suffix or any input suffix
if options[:phylip_suffix]
  if options[:phylip_suffix] == options[:output_suffix]
    abort "\n  ERROR: phylip suffix and output suffix are identical ('#{options[:phylip_suffix]}').\n" \
          "  Use -p and -x to set different suffixes.\n"
  end
  colliding = [options[:input_suffix], (options[:aligned_suffix] if options[:add_mode])].compact
  colliding.each do |s|
    if options[:phylip_suffix] == s
      abort "\n  ERROR: phylip suffix ('#{options[:phylip_suffix]}') matches an input suffix.\n" \
            "  This would overwrite source files. Use a different suffix for -p.\n"
    end
  end
end

# Check for MAFFT — distinguish "not found" from other failures
mafft_check = system('which mafft > /dev/null 2>&1')
if mafft_check.nil?
  abort "\n  ERROR: could not execute 'which' to locate mafft.\n"
elsif !mafft_check
  abort "\n  ERROR: 'mafft' executable not found in PATH.\n"
end

mafft_version = `mafft --version 2>&1`.lines.first&.strip || 'unknown'

# ── Build effective MAFFT options ─────────────────────────────────────────────

mafft_opts = options[:mafft_opts].dup
mafft_opts += " --thread #{options[:threads]}" if options[:threads]

# In add mode, --keeplength is positional: it must appear between --add <unaligned>
# and <aligned> in the command, not in the general opts string.  Strip it here so
# that a user passing it via -m cannot accidentally misplace or duplicate it.
mafft_opts = mafft_opts.gsub(/\s*--keeplength\b/, '').strip if options[:add_mode]

# ── File discovery ────────────────────────────────────────────────────────────

if options[:add_mode]
  # Scan for already-aligned files; match each with its unaligned counterpart by stem
  aligned_files = Dir.glob(File.join(options[:input_dir], "*#{options[:aligned_suffix]}")).sort

  if aligned_files.empty?
    abort "\n  ERROR: no files matching *#{options[:aligned_suffix]} found in '#{options[:input_dir]}'.\n"
  end

  pairs         = []
  missing_pairs = []

  aligned_files.each do |aligned_file|
    stem          = File.basename(aligned_file, options[:aligned_suffix])
    unaligned_file = File.join(options[:input_dir], "#{stem}#{options[:input_suffix]}")

    if File.exist?(unaligned_file)
      pairs << { stem: stem, aligned: aligned_file, unaligned: unaligned_file }
    else
      missing_pairs << { aligned: aligned_file, expected: unaligned_file }
    end
  end

  if missing_pairs.any?
    puts "  WARNING: #{missing_pairs.size} aligned file(s) have no matching unaligned file and will be skipped:"
    missing_pairs.each { |mp| puts "    • #{File.basename(mp[:aligned])}  (expected #{File.basename(mp[:expected])})" }
    puts
  end

  if pairs.empty?
    abort "\n  ERROR: no matched pairs found. Check your suffixes.\n"
  end

else
  files = Dir.glob(File.join(options[:input_dir], "*#{options[:input_suffix]}")).sort

  abort "\n  ERROR: no files matching *#{options[:input_suffix]} found in '#{options[:input_dir]}'.\n" \
    if files.empty?
end

# ── Summary ───────────────────────────────────────────────────────────────────

work_items = options[:add_mode] ? pairs : files.map { |f| { stem: File.basename(f, options[:input_suffix]), file: f } }

skipped_count = 0
unless options[:overwrite]
  skipped_count = work_items.count do |item|
    stem        = options[:add_mode] ? item[:stem] : File.basename(item[:file], options[:input_suffix])
    output_file = File.join(options[:output_dir], "#{stem}#{options[:output_suffix]}")
    File.exist?(output_file)
  end
end

puts "  Mode         : #{options[:add_mode] ? 'add (--add --keeplength)' : 'standard'}"
puts "  Input dir    : #{options[:input_dir]}"
if options[:add_mode]
  puts "  Aligned suf  : #{options[:aligned_suffix]}"
  puts "  Unaligned suf: #{options[:input_suffix]}"
  puts "  Matched pairs: #{pairs.size}"
else
  puts "  Input suffix : #{options[:input_suffix]}"
  puts "  Files found  : #{files.size}"
end
puts "  Output dir   : #{options[:output_dir]}"
puts "  Output suffix: #{options[:output_suffix]}"
puts "  MAFFT opts   : #{mafft_opts}"
puts "  --keeplength : always enforced in add mode" if options[:add_mode]
puts "  Threads      : #{options[:threads] || '(MAFFT default)'}"
puts "  Overwrite    : #{options[:overwrite]}"
puts "  Phylip suffix: #{options[:phylip_suffix]}" if options[:phylip_suffix]
puts "  MSA quality  : #{options[:msa_quality]}"   if options[:msa_quality]
puts "  Dry run      : yes — no files will be written" if options[:dry_run]
puts "  Will skip    : #{skipped_count} (output already exists)" unless options[:overwrite] || skipped_count.zero?
puts "  Log file     : #{options[:log_file_path]}" if options[:log_file_path]
puts

# ── Open log file and write header ────────────────────────────────────────────

run_started_at = Time.now

log_io = nil
if options[:log_file_path] && !options[:dry_run]
  FileUtils.mkdir_p(File.dirname(File.expand_path(options[:log_file_path])))
  log_io = File.open(options[:log_file_path], 'w')

  log_io.puts '# MAFFT Batch Aligner — run log'
  log_io.puts "# Script version : #{VERSION}"
  log_io.puts "# MAFFT version  : #{mafft_version}"
  log_io.puts "# Started        : #{run_started_at.strftime('%Y-%m-%d %H:%M:%S %Z')}"
  log_io.puts "# Command        : #{File.basename($PROGRAM_NAME)} #{original_argv.join(' ')}"
  log_io.puts '#'
  log_io.puts "# Mode           : #{options[:add_mode] ? 'add (--add --keeplength)' : 'standard'}"
  log_io.puts "# Input dir      : #{options[:input_dir]}"
  if options[:add_mode]
    log_io.puts "# Aligned suf    : #{options[:aligned_suffix]}"
    log_io.puts "# Unaligned suf  : #{options[:input_suffix]}"
    log_io.puts "# Matched pairs  : #{pairs.size}"
  else
    log_io.puts "# Input suffix   : #{options[:input_suffix]}"
    log_io.puts "# Files found    : #{files.size}"
  end
  log_io.puts "# Output dir     : #{options[:output_dir]}"
  log_io.puts "# Output suffix  : #{options[:output_suffix]}"
  log_io.puts "# MAFFT opts     : #{mafft_opts}"
  log_io.puts "# --keeplength   : always enforced in add mode" if options[:add_mode]
  log_io.puts "# Threads        : #{options[:threads] || '(MAFFT default)'}"
  log_io.puts "# Overwrite      : #{options[:overwrite]}"
  log_io.puts "# Phylip suffix  : #{options[:phylip_suffix]}" if options[:phylip_suffix]
  log_io.puts "# MSA quality    : #{options[:msa_quality]}"   if options[:msa_quality]
  log_io.puts '#'
  log_io.puts "# #{%w[stem status elapsed_s input output reason].join("\t")}"
  log_io.flush
end

# ── Open MSA quality file ─────────────────────────────────────────────────────

MSA_QUALITY_COLS = %w[
  stem ntax nchar
  gap_pct mean_seq_gap_pct
  max_gap_pct max_gap_seq
  min_gap_pct min_gap_seq
  conserved_cols all_gap_cols
  cols_gt50pct_gaps cols_gt50pct_gaps_pct
  cols_gt90pct_gaps cols_gt90pct_gaps_pct
].freeze

qual_io = nil
if options[:msa_quality] && !options[:dry_run]
  FileUtils.mkdir_p(File.dirname(File.expand_path(options[:msa_quality])))
  qual_io = File.open(options[:msa_quality], 'w')
  qual_io.puts MSA_QUALITY_COLS.join("\t")
  qual_io.flush
end

# ── Run MAFFT ─────────────────────────────────────────────────────────────────

bar     = options[:dry_run] ? nil : ProgressBar.new(work_items.size)
errors  = []
skipped = []

work_items.each do |item|
  if options[:add_mode]
    stem          = item[:stem]
    aligned_file  = item[:aligned]
    unaligned_file = item[:unaligned]
    label         = File.basename(aligned_file)
  else
    stem          = File.basename(item[:file], options[:input_suffix])
    input_file    = item[:file]
    label         = File.basename(input_file)
  end

  output_file = File.join(options[:output_dir], "#{stem}#{options[:output_suffix]}")
  log_file    = File.join(options[:output_dir], "#{stem}.mafft.log")
  input_desc  = options[:add_mode] ? "#{aligned_file} + #{unaligned_file}" : input_file

  # Skip existing outputs unless --overwrite
  unless options[:overwrite]
    if File.exist?(output_file)
      skipped << label
      if options[:dry_run]
        puts "  [SKIP] #{label}\n         output already exists: #{output_file}"
      else
        log_io&.puts [stem, 'skipped', 0, input_desc, output_file, 'output already exists'].join("\t")
        log_io&.flush
        bar.increment(label)
      end
      next
    end
  end

  # ── FASTA content validation ─────────────────────────────────────────────────

  issues = if options[:add_mode]
    v  = validate_fasta(aligned_file, 'aligned')
    v += check_equal_lengths(aligned_file, 'aligned') if v.empty?
    if v.empty?
      v += validate_fasta(unaligned_file, 'unaligned')
      v += check_id_overlap(aligned_file, unaligned_file) if v.empty?
    end
    v
  else
    validate_fasta(input_file, 'input')
  end

  unless issues.empty?
    reason = issues.join('; ')
    errors << { file: label, reason: reason }
    if options[:dry_run]
      puts "  [FAIL] #{label}\n         #{reason}"
    else
      log_io&.puts [stem, 'failed', 0, input_desc, '-', reason].join("\t")
      log_io&.flush
      bar.increment(label)
    end
    next
  end

  # ── Build MAFFT command ───────────────────────────────────────────────────────

  cmd = if options[:add_mode]
    "mafft #{mafft_opts} --add #{Shellwords.escape(unaligned_file)} --keeplength " \
    "#{Shellwords.escape(aligned_file)} " \
    "> #{Shellwords.escape(output_file)} " \
    "2> #{Shellwords.escape(log_file)}"
  else
    "mafft #{mafft_opts} #{Shellwords.escape(input_file)} " \
    "> #{Shellwords.escape(output_file)} " \
    "2> #{Shellwords.escape(log_file)}"
  end

  if options[:dry_run]
    puts "  [RUN]  #{label}\n         #{cmd}"
    if options[:phylip_suffix]
      puts "  [PHY]  #{label}\n         #{File.join(options[:output_dir], "#{stem}#{options[:phylip_suffix]}")}"
    end
    puts "  [MSQ]  #{label}\n         metrics → #{options[:msa_quality]}" if options[:msa_quality]
    next
  end

  # ── Execute ───────────────────────────────────────────────────────────────────

  file_start = Time.now
  result     = system(cmd)
  file_secs  = (Time.now - file_start).round(2)

  if result.nil?
    reason = 'mafft could not be executed (command not found?)'
    errors << { file: label, reason: reason }
    FileUtils.rm_f(output_file)
    log_io&.puts [stem, 'failed', file_secs, input_desc, '-', reason].join("\t")
  elsif !result
    exit_code = $CHILD_STATUS.exitstatus
    reason    = "mafft exited with status #{exit_code} (see #{File.basename(log_file)})"
    errors << { file: label, reason: reason }
    FileUtils.rm_f(output_file)
    log_io&.puts [stem, 'failed', file_secs, input_desc, '-', reason].join("\t")
  else
    log_io&.puts [stem, 'succeeded', file_secs, input_desc, output_file, ''].join("\t")

    # Parse the output once; reused for PHYLIP conversion and quality metrics
    aligned_records = (options[:phylip_suffix] || qual_io) ? parse_fasta(output_file) : []

    if options[:phylip_suffix]
      phylip_file  = File.join(options[:output_dir], "#{stem}#{options[:phylip_suffix]}")
      phy_warnings = write_phylip(aligned_records, phylip_file)
      phy_warnings.each do |w|
        log_io&.puts [stem, 'warning', '-', '-', phylip_file, w].join("\t")
        $stderr.puts "  WARNING (#{label}): #{w.strip}"
      end
    end

    if qual_io
      m = compute_msa_metrics(aligned_records)
      if m
        row = MSA_QUALITY_COLS.drop(1).map { |k| m[k.to_sym] }   # drop 'stem'; prepend below
        qual_io.puts ([stem] + row).join("\t")
        qual_io.flush
      end
    end
  end
  log_io&.flush

  bar.increment(label)
end

bar&.finish

# ── Final report ──────────────────────────────────────────────────────────────

elapsed_str = options[:dry_run] ? 'n/a' : format_elapsed(bar.elapsed)
succeeded   = work_items.size - errors.size - skipped.size

puts "  Results"
puts "  ───────────────────────────────────────────"
puts "  Total        : #{work_items.size}"
puts "  Succeeded    : #{succeeded}"
puts "  Skipped      : #{skipped.size} (already existed)" unless skipped.empty?
puts "  Failed       : #{errors.size}"
puts "  Elapsed time : #{elapsed_str}"

if errors.any?
  puts
  puts "  Failed files:"
  errors.each { |e| puts "    • #{e[:file]}\n      reason: #{e[:reason]}" }
end

if log_io
  log_io.puts '#'
  log_io.puts '# Summary'
  log_io.puts "# Total     : #{work_items.size}"
  log_io.puts "# Succeeded : #{succeeded}"
  log_io.puts "# Skipped   : #{skipped.size}" unless skipped.empty?
  log_io.puts "# Failed    : #{errors.size}"
  log_io.puts "# Elapsed   : #{elapsed_str}"
  log_io.puts "# Finished  : #{Time.now.strftime('%Y-%m-%d %H:%M:%S %Z')}"
  log_io.close
  puts "  Log written  : #{options[:log_file_path]}"
end

if qual_io
  qual_io.close
  puts "  Quality file : #{options[:msa_quality]}"
end

puts

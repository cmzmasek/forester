#
# To change this template, choose Tools | Templates
# and open the template in the editor.


class HmmscanParser

  def initialize file
    @file = file
  end

  def parse
    results = []
    File.open( @file ).each do | line |
      line.strip!
      if !line.empty? && line[ 0 ] != "#"
        results << parse_line( line )
      end
    end
    results
  end

  def parse_line( line )
    #         tn      acc     tlen    query   acc     qlen    Evalue  score   bias    #       of      c-E     i-E     score   bias    hf      ht      af      at      ef      et      acc     desc
    #         1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21      22      23
    if line =~ /^(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(.*)/

      r = HmmscanResult.new

      r.model      = $1
      r.tlen       = $3.to_i
      r.query      = $4
      r.qlen       = $6.to_i

      r.fs_e_value = $7.to_f
      r.fs_score   = $8.to_f

      r.number     = $10.to_i
      r.out_of     = $11.to_i

      r.c_e_value  = $12.to_f
      r.i_e_value  = $13.to_f
      r.score      = $14.to_f

      r.hmm_from   = $16.to_i
      r.hmm_to     = $17.to_i

      r.ali_from   = $18.to_i
      r.ali_to     = $19.to_i

      r.env_from   = $20.to_i
      r.env_to     = $21.to_i

      if r.number > r.out_of || r.hmm_from > r.hmm_to || r.ali_from > r.ali_to || r.env_from > r.env_to
        raise IOError, "illogical format: " + line
      end

      r
    else
      raise IOError, "illegal format: " + line
    end
  end
end

class HmmscanResult

  attr_accessor :model
  attr_accessor :tlen
  attr_accessor :query
  attr_accessor :qlen
  attr_accessor :fs_e_value
  attr_accessor :fs_score
  attr_accessor :number
  attr_accessor :out_of
  attr_accessor :c_e_value
  attr_accessor :i_e_value
  attr_accessor :score
  attr_accessor :hmm_from
  attr_accessor :hmm_to
  attr_accessor :ali_from
  attr_accessor :ali_to
  attr_accessor :env_from
  attr_accessor :env_to

end

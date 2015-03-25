require 'stringio'

module BioMummer

  class Alignment
    attr_accessor :refname, :qryname, :refstart, :refstop, :qrystart, :qrystop, :strand, :distances

    def initialize(refname, qryname, refstart, refstop, qrystart, qrystop, strand, distances)
      @refname = refname
      @qryname = qryname
      @refstart = refstart
      @refstop = refstop
      @qrystart = qrystart
      @qrystop = qrystop
      @strand = strand
      @distances = distances
    end

    def deltas
      ds = []
      a = @distances.each_with_object([0]) do |d, arr|
        state = arr.last
        if d > 0
          ds += Array.new(d - 1, state)
          ds.push(nil)
          arr << state - 1
        else
          ds += Array.new(d * -1 - 1, state)
          arr << state + 1
        end
      end
      ds << a.last
      return ds
    end

    def ref_to_query(ref_position)
      position_in_alignment = ref_position - refstart + 1
      qrypos = nil
      if position_in_alignment >= deltas.length
        qrypos = @qrystart - 1 + position_in_alignment + deltas.last
      elsif deltas[position_in_alignment-1]
        qrypos = @qrystart - 1 + position_in_alignment + deltas[position_in_alignment-1]
      end
      !@strand && qrypos ? @qrystart + @qrystop - qrypos : qrypos
    end
  end
  
  class DeltaFile
    NUCMER = 1
    PROMER = 2

    attr_accessor :reference_filename, :alternate_filename
    attr_reader :alignments
    
    def initialize(io)
      raise EncodingError, "Not delta format" unless format_ok?(io)
      @alignments = parse(io.read)
    end

    def self.open(filename)
      io = File.open(filename)
      self.new(super(io))
    end

    def parse(string)
      string.split("\n").slice_before(/^>/).flat_map do |block|
        refname, qryname = block.shift.match(/>(.*) (.*) \d+ \d+/).captures
        block.slice_before(/\d+ \d+ \d+ \d+/).map do |alignment|
          refstart, refstop, qrystart, qrystop = alignment
            .shift
            .match(/(\d+) (\d+) (\d+) (\d+) /)
            .captures
            .map{ |c| c.to_i }
          alignment.pop
          Alignment.new(refname,
                        qryname,
                        refstart,
                        refstop,
                        [qrystart, qrystop].min,
                        [qrystart, qrystop].max,
                        qrystart < qrystop,
                        alignment.map{ |i| i.to_i })
        end
      end
    end

    def format_ok?(io)
      line = io.gets.chomp
      unless line.match(/^(\/.*) (\/.*)$/)
        return false
      else
        reference_filename = $1
        alternate_filename = $2
      end

      case io.gets.chomp
      when /NUCMER/
        @format = NUCMER
        return true
      when /PROMER/
        @format = PROMER
        return true
      else
        return false
      end
    end

    def transpose_region(refname, startpos, endpos)
      a = alignments.find do |a|
        a.refname == refname && a.refstart <= startpos && a.refstop >= endpos
      end
      if a
        qryname = a.qryname
        qrystart = a.ref_to_query(startpos)
        qrystop = a.ref_to_query(endpos)
        if qrystart.nil? || qrystop.nil?
          return nil
        else
          return [qryname, qrystart, qrystop, a.strand]
        end
      else
        return nil
      end
    end

    def self.open(filename)
      return new(File.open(filename, 'r'))
    end
  end
end

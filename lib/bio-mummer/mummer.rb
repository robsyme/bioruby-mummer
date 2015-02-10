
module BioMummer

  class Alignment
    attr_accessor :refname, :queryname, :refstart, :refstop, :querystart, :querystop, :strand, :distances

    def initialize(refname, queryname, refstart, refstop, querystart, querystop, strand, distances)
      @refname = refname
      @queryname = queryname
      @refstart = refstart
      @refstop = refstop
      @querystart = querystart
      @querystop = querystop
      @strand = strand
      @distances = distances
      @deltas = []
    end

    def deltas
      if @deltas == []
        a = @distances.each_with_object([0]) do |d, arr|
          state = arr.last
          if d > 0
            @deltas += Array.new(d - 1, state)
            @deltas.push(nil)
            arr << state - 1
          else
            @deltas += Array.new(d * -1 - 1, state)
            arr << state + 1
          end
        end
        @deltas << a.last
      end
      return @deltas
    end

    def ref_to_query(ref_position)
      position_in_alignment = ref_position - refstart + 1
      querypos = nil
      if position_in_alignment >= deltas.length
        querypos = @querystart - 1 + position_in_alignment + deltas.last
      elsif deltas[position_in_alignment-1]
        querypos = @querystart - 1 + position_in_alignment + deltas[position_in_alignment-1]
      end
      !@strand && querypos ? @querystart + @querystop - querypos : querypos
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

    def parse(string)
      string.split("\n").slice_before(/^>/).flat_map do |block|
        refname, queryname = block.shift.match(/>(.*) (.*) \d+ \d+/).captures
        block.slice_before(/\d+ \d+ \d+ \d+/).map do |alignment|
          refstart, refstop, querystart, querystop = alignment
            .shift
            .match(/(\d+) (\d+) (\d+) (\d+) /)
            .captures
            .map{ |c| c.to_i }
          alignment.pop
          Alignment.new(refname,
                        queryname,
                        refstart,
                        refstop,
                        [querystart, querystop].min,
                        [querystart, querystop].max,
                        querystart < querystop,
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
        queryname = a.queryname
        querystart = a.ref_to_query(startpos)
        querystop = a.ref_to_query(endpos)
        if querystart.nil? || querystop.nil?
          return nil
        else
          return [queryname, querystart, querystop, a.strand]
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

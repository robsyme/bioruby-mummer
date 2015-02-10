require 'helper'

class TestBioMummer < MiniTest::Test

  def setup
    @report = BioMummer::DeltaFile.open("test/data/out.delta")
  end
  
  should "open and close a delta file" do
    assert_kind_of BioMummer::DeltaFile, @report
  end

  should "have enumerable alignments" do
    assert_equal "reference", @report.alignments.first.refname
    assert_kind_of Fixnum, @report.alignments.first.refstart
    assert_equal 1, @report.alignments.first.refstart
    assert_equal 2435, @report.alignments.first.refstop
    assert_equal 1, @report.alignments.first.querystart
    assert_equal 2435, @report.alignments.first.querystop
  end

  should "do basic coordinate transforms on the positive strand" do
    a = BioMummer::Alignment.new('ref', 'query', 101, 150, 201, 248, true, [5,-2,3,1,1,-4])
    assert_equal 201, a.ref_to_query(101)
    assert_equal 202, a.ref_to_query(102)
    assert_equal 203, a.ref_to_query(103)
    assert_equal 204, a.ref_to_query(104)
    assert_equal nil, a.ref_to_query(105)
    assert_equal 205, a.ref_to_query(106)
    assert_equal 207, a.ref_to_query(107)
    assert_equal 208, a.ref_to_query(108)
    assert_equal nil, a.ref_to_query(109)
    assert_equal nil, a.ref_to_query(110)
    assert_equal nil, a.ref_to_query(111)
    assert_equal 209, a.ref_to_query(112)
    assert_equal 217, a.ref_to_query(119)
    assert_equal 248, a.ref_to_query(150)
  end

  should "do basic coordinate transforms on the negative strand" do
    a = BioMummer::Alignment.new('ref', 'query', 101, 118, 201, 216, false, [5,-2,3,1,1,-4])
    assert_equal 216, a.ref_to_query(101)
    assert_equal 201, a.ref_to_query(118)
  end
end


classdef SourcePrimary 
  properties
    rumi
    rice
    wast
    gas
    coal
    rco
    otherff
    bb
    wet
    trmt
    anim
    geo
  end

  methods(Static)
    function inout = clearMemory(inout)
      for name = ["rumi","rice","wast","gas","coal","rco","otherff","wet","trmt","anim"]
        inout.(name) = [];
      end
    end
  end
end

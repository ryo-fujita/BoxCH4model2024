% Handles calculation for global source fractions
% Ryo Fujtia 2024
classdef SourceFraction
  properties
    tspan 
    bio  
    ff
    bb
    anth_bio 
    natr_bio 
    anth_ff
    geo
  end

  methods
    function obj = SourceFraction(varargin) 
      if nargin ~= 0
        model = varargin{1};
        obj(size(model.para)) = obj; 

        if nargin == 2 && ~isempty(varargin{2})
          nameList = varargin{2};
          emsCH4   = Emission().getCH4only(model,[nameList,"tot"]); 
        else
          nameList = ["bio","ff","bb","anth_bio","natr_bio","anth_ff","geo"];
          emsCH4   = Emission().getCH4only(model,[nameList,"tot"]); 
        end

        for iRun = 1:size(obj,1)
          for iNode = 1:size(obj,2)
            obj(iRun,iNode).tspan = emsCH4(iRun,iNode).tspan;
            for name = nameList
              obj(iRun,iNode).(name) = emsCH4(iRun,iNode).(name)./emsCH4(iRun,iNode).tot .* 100;
            end
          end
        end
      end
    end
  end

end
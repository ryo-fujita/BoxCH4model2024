% Clear memory
% Ryo Fujita 2024

classdef ClearMemory
  methods (Static)
    function y = getFieldList(obj,varargin)
      if ~isempty(varargin), nameList = varargin{1}; else, nameList = []; end 

      nameListAll = string(fieldnames(obj)).';
      indexClear  = ~ismember(nameListAll,["tspan" nameList]); 
      y = nameListAll(indexClear);
    end

%%
    function obj = Select(obj,varargin)
      if ~isempty(varargin)
        clearNameList = varargin{1};
      else
      switch class(obj)
        case 'Atmosphere'
          clearNameList = ["C13","CHD","C14"];
        case 'Emission'
          clearNameList = ["C13","CHD","C14"]; 
        case 'Loss'
          clearNameList = ["C13","CHD","C14"];
      end
      end

      for iRun = 1:size(obj,1)
      for iNode= 1:size(obj,2)
        for nameID = clearNameList
          obj(iRun,iNode).(nameID) = [];
        end
      end  
      end
    end

%% 
    function inputObj = SelectStruct(inputObj,propList,structList)
      for iNode = 1:length(inputObj)
        for pname = propList
          if isempty(inputObj(iNode).(pname)), continue; end 
          inputObj(iNode).(pname) = rmfield(inputObj(iNode).(pname), structList);
        end
      end
    end

  end
end


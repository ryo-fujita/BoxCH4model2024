% Control title in figure
% Ryo Fujita 2024

classdef Title < Plot
  properties
    nameClass 
    filename
    strFigTop
  end
%%
  properties (SetAccess = protected)
    calcType       
    dataTypeList   
    statsNameList  
  end

  methods
    function obj = Title(workdirList,model,plotType,categoryList,...
      periodList,lgdNameList,calcType,statsNameList,nameClass,varargin)

      obj = obj@Plot(workdirList,model,plotType,categoryList); %Explicitly call the parent constructor first
      obj.lgdNameList   = lgdNameList;
      obj.periodList    = periodList;
      obj.calcType      = calcType;
      obj.statsNameList = statsNameList;
      obj.nameClass     = nameClass;
      
      if size(obj.categoryList,1)  > 1, obj.categoryList  = obj.categoryList.'; end
      if size(obj.lgdNameList,1)   > 1, obj.lgdNameList   = obj.lgdNameList.'; end
      if size(obj.statsNameList,1) > 1, obj.statsNameList = obj.statsNameList.'; end
      
      switch obj.nameClass 
        case 'CompareScenario'
          obj.filename = obj.filenameComp;
        case 'Gplotmatrix' 
          plotTypeY     = varargin{1};
          categoryListY = varargin{2};
          periodListY   = varargin{3};
          calcTypeY     = varargin{4};
          obj.filename  = obj.filenameGplot(plotTypeY,categoryListY,periodListY,calcTypeY);
        otherwise
          obj.filename  = obj.filenameDef; 
      end
      
      if length(char(obj.filename)) > 240 
        titleChar    = char(obj.filename);
        obj.filename = titleChar(1:240);
      end
    end  

%% 
    function y = filenameDef(obj)
      plotType = cutLength(obj.plotType,9);
      period   = mergeArrayToOne(obj.periodList);
      stats    = mergeArrayToOne(obj.statsNameList);
      category = mergeArrayToOne(obj.categoryList);
      lgd      = mergeArrayToOne(obj.lgdNameList);

      y = strcat(obj.nameClass,plotType,obj.calcType,period,stats,category,lgd); 
    end

%% 
    function y = filenameGplot(obj,plotTypeY,categoryListY,periodListY,calcTypeY)
      plotTypeX = cutLength(obj.plotType,9);
      periodX   = mergeArrayToOne(obj.periodList);
      categoryX = mergeArrayToOne(obj.categoryList);
      lgd       = mergeArrayToOne(obj.lgdNameList);
      
      if isempty(plotTypeY)
        y = strcat(obj.nameClass,plotTypeX,obj.calcType,periodX,categoryX,lgd); 
      else
        plotTypeY = cutLength(plotTypeY,9); 
        periodY   = mergeArrayToOne(periodListY);
        categoryY = mergeArrayToOne(categoryListY);
        if strcmp(plotTypeY,plotTypeX),    plotTypeY = ""; end
        if strcmp(periodY,periodX),        periodY   = ""; end
        if strcmp(categoryY,categoryX),    categoryY = ""; end
        if strcmp(calcTypeY,obj.calcType), calcTypeY = ""; end
        if length(obj.workdirList) > 1, flag = "Mean"; else, flag = ""; end
        y = strcat(obj.nameClass,flag,plotTypeX,obj.calcType,periodX, ...
          categoryX,"_",plotTypeY,calcTypeY,periodY,categoryY,lgd); 
      end
    end

%%
    function y = filenameComp(obj)
      strPrdComb = ""; 
      for iPrd = length(obj.periodList)
        for iRow = 1:size(obj.periodList{iPrd},1)
          strPrdComb = strcat(strPrdComb,"_",num2str(obj.periodList{iPrd}(iRow,1)), ...
            "-",num2str(obj.periodList{iPrd}(iRow,2))); 
        end
      end

      strInvComb = ""; for strInv = obj.lgdNameList,   strInvComb = strcat(strInvComb,"_",strInv); end
      strPrcComb = ""; for strPrc = obj.statsNameList, strPrcComb = strcat(strPrcComb,strPrc); end
      strCatComb = ""; for strCat = obj.categoryList,  strCatComb = strcat(strCatComb,"_",strCat); end
      y = strcat(obj.plotType,obj.calcType,strPrdComb,strPrcComb,strCatComb,strInvComb); 
    end

%% Not used but needed as an abstrct class for Plot.m
    function postDo(obj)
    end  
    function plotDataCategory(obj)
    end  
    function setupXYLabelPrep(obj)
    end  
  end
end

%% 
function y = cutLength(str,numLim)
  len = length(char(str)); if len > numLim, len = numLim; end
  y = extractBefore(str,len+1);
end

%%
function y = mergeArrayToOne(array)
  y = ""; 
  switch class(array)
  case 'cell' %periodList
    if size(array,1) ~= 1 && size(array,2) ~= 2
      error('cell size(array,1) ~= 1! not proper for the title of Prctile'); 
    end

    for i = length(array)
      for iRow = 1:size(array{i},1)
        array{i}(iRow,1) = floor(array{i}(iRow,1)); 
        array{i}(iRow,2) = floor(array{i}(iRow,2));

        if isequal(array{i}(iRow,1),array{i}(iRow,2))
          y = strcat(y,"_",num2str(array{i}(iRow,1)));
        else
          y = strcat(y,"_",num2str(array{i}(iRow,1)),"-",num2str(array{i}(iRow,2)));
        end
      end
    end
  otherwise
    for str = array, y = strcat(y,"_",str); end
  end
end
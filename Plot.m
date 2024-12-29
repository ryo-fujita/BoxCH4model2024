% An abstract superclass that aggregates the properties and 
% functions for figure processing. Note that not all the children 
% functions are openly available. Please contact the author if you 
% would like to use unavailable plot functions. 
% Ryo Fujita 2024

classdef Plot < handle
%%
  properties (SetAccess = private)
    workdirList string = "."; 
  end

%%
  properties
    model Model
    modStyle ModifyPlotStyle
    plotObj       

    workdir       string
    lgdNameList   string
    plotType      string
    categoryList  string
    periodList    cell 
    subPlotLayout double
    inventory     string
    lagTermList   double 

    iDir            double = 1
    currentCategory string
    period          double
    currentDataType string

    flagScenarioMean string = 'n'; %'y' or 'n'
    setSaveFig  string = 'yes'; 
    setMean     string ; %[] or "aveHist" or "aveWT"

    plotClass string

    optFig struct = struct('savePNG',true,'saveFIG',false,'close',false);
  end

%%
  properties (Abstract, SetAccess = protected)
    calcType      string
    dataTypeList  string 
    statsNameList string
  end

%%
  methods 
    function obj = Plot(workdirList,model,plotType,categoryList,varargin)
      if isa(obj,'CompareScenario')
        return; 
      end

      if nargin > 0
        obj.workdirList  = workdirList;
        obj.model     = model;
        obj.plotType     = plotType;
        obj.categoryList = categoryList;
        obj.addProp(varargin{:});
      end

      obj.checkEmpty("model");
      obj.checkEmpty("plotType");
      obj.checkEmpty("categoryList");
      obj.checkEmpty("plotObj"); 

      obj.plotClass = class(obj);
    end

%%
    function addProp(obj,varargin)
      if length(varargin) < 2, return; end
      for k=1:2:length(varargin)
        obj.(varargin{k}) = varargin{k+1};
      end
    end

%%
    function h = plot(obj,varargin)
      obj.addProp(varargin{1:length(varargin)});
      obj.checkEmpty("periodList");  
      obj.checkEmpty("lgdNameList");
      obj.checkEmpty("calcType"); 
      obj.checkEmpty("dataTypeList"); 
      obj.checkEmpty("statsNameList"); 
      obj.checkEmpty("subPlotLayout"); 
      
      if isa(obj,'Histogram'), obj.updateHistCommonLim; end

      figure('Name',strcat(class(obj),obj.plotType,obj.calcType)); 
      obj.modStyle  = ModifyPlotStyle(obj.model);

      for i = 1:length(obj.workdirList) 
        obj.iDir = i;
        obj.workdir   = obj.workdirList(obj.iDir);
        obj.model  = obj.model.setup(obj.workdir);
        obj.inventory = obj.lgdNameList(obj.iDir);
        obj.plotData;
      end
      obj.postDo; 
      obj.saveFigure; h = gcf;
      obj.clearProp(["subPlotLayout"]);
    end

%%
    function obj = plotData(obj)
      for iCategory = 1:length(obj.categoryList)
        obj.currentCategory = obj.categoryList(iCategory);
        obj.subplotFront(iCategory);
        obj.setupXYLabelPrep;       
        obj.plotDataCategory; 
        pause(0.2); %required here to stabilize the graphic system
      end
    end

%% 
    function subplotFront(obj,iCategory)
      if length(obj.categoryList) == 1, return; end
      subplot(obj.subPlotLayout(1),obj.subPlotLayout(2),iCategory);
    end

%% 
     function saveFigure(obj)  
      if ~strcmp(obj.setSaveFig,'yes'), return; end 
      switch class(obj)
        case 'Gplotmatrix'
          title = Title(obj.workdirList,obj.model,obj.plotType,...
            obj.categoryList,obj.periodList,obj.lgdNameList,obj.calcType, ...
            obj.statsNameList,class(obj),obj.plotTypeY,obj.categoryListY, ...
            obj.periodListY,obj.calcTypeY); 
        otherwise
          title = Title(obj.workdirList,obj.model,obj.plotType, ...
            obj.categoryList,obj.periodList,obj.lgdNameList,obj.calcType, ...
            obj.statsNameList,class(obj)); 
      end

      if obj.optFig.savePNG
        saveas(gcf,strcat(obj.workdir,"/",title.filename,'.png')); 
      end 

      if obj.optFig.saveFIG
        saveas(gcf,strcat(obj.workdir,"/",title.filename,'.fig')); 
      end 
      
      if obj.optFig.close
        close(gcf); 
      end 
    end

%%
    function h = plotMean(obj,varargin)
      obj.mean;
      if nargin >= 3
        for k=1:2:length(varargin)
          obj.(varargin{k}) = varargin{k+1};
        end
      end
      h = obj.plot;
    end

%%
    function obj = mean(obj)
      if size(obj.plotObj,2) < 2, error('size(obj.plotObj,2) < 2!'); end  
      obj.setMean = "aveHist"; 
      obj.plotObj = obj.statsHistMean(obj.plotObj,obj.categoryList,obj.periodList);

      obj.workdirList = obj.workdirList(1);
      obj.lgdNameList = obj.setMean; 
    end

%% 
    function y = statsHistMean(obj,plotObj,nameList,periodList)
      obj.checkErrorHistMean(plotObj,nameList);
      y.classType = class(plotObj(1)); 
      y.tspan     = plotObj(1).tspan;

      if ~strcmp(obj.calcType,'Timeseries')
        for i = 1:size(plotObj,2)
          plotObj(i) = MyCalc.convertCalcType(obj.calcType,plotObj(i),nameList,periodList);
        end
      end

      for name = nameList
        y.(name) = obj.calcStatsHistMean(name,obj.plotObj3D(plotObj,name));
      end
    end

%% 
    function y = calcStatsHistMean(obj,name,obj3D)
      for i = 1:size(obj3D,1) %i: index for time
        y.NumBins        = 500;
        y.BinLimits(i,:) = MyCalc.getBinLimitsList(...
          obj.plotType,obj.calcType,name,obj.model,squeeze(obj3D(i,:,:))).(name); 
        y.BinEdges(i,:)  = MyCalc.binEdges(y.BinLimits(i,:),y.NumBins);
        y.Values(i,:)    = MyCalc.meanPropability(squeeze(obj3D(i,:,:)),y.BinEdges(i,:));

        y.cdf(i,:) = [0 cumsum(y.Values(i,:))]; 
        y.ave(i)   = sum(movmean(y.BinEdges(i,:),2,'Endpoints','discard').*y.Values(i,:));
        [~,idx1]   = min(abs(y.cdf(i,:)-0.5));
        y.med(i)   = y.BinEdges(i,idx1);
        y.cdf50(i) = y.cdf(i,idx1);

        [~,idx1] = min(abs(y.cdf(i,:)-0.16));
        [~,idx2] = min(abs(y.cdf(i,:)-0.84));
        y.prc68(i,:) = [y.BinEdges(i,idx1) y.BinEdges(i,idx2)];
        y.cdf68(i,:) = [y.cdf(i,idx1) y.cdf(i,idx2)];

        [~,idx1] = min(abs(y.cdf(i,:)-0.025));
        [~,idx2] = min(abs(y.cdf(i,:)-0.975));
        y.prc95(i,:) = [y.BinEdges(i,idx1) y.BinEdges(i,idx2)];
        y.cdf95(i,:) = [y.cdf(i,idx1) y.cdf(i,idx2)];
      end

      y.NumCase_org = ones(1,size(obj3D,3)).*size(obj3D,2); 

      % for save the strage of memory
      if strcmp(obj.calcType,'Timeseries') 
        clearNameList = ClearMemory.getFieldList(y,["classType","ave","med","prc68","prc95"]);
        y = ClearMemory.Select(y,clearNameList);
      end

      y.plotStyleType = 'Histogram';
    end

%% Name of legend after averaging multiple inventories' results
    function setMeanInfo(obj,val)
      switch val
      case {'aveHist','aveWT'}
        obj.setMean = val;
      otherwise
        error(strcat(val,' is not defined as obj.setMean!'))
      end
    end

%%
    function checkEmpty(obj,propName)
      if sum(ismember(string(fieldnames(obj)),propName),'all') == 0
        error(strcat('Undefined input propertiy name: ',propName))
      elseif isempty(obj.(propName))
        switch propName
        case 'model'
          obj.(propName) = Model(obj.workdirList(1),"Posterior",9); 
        case 'lgdNameList'
          obj.(propName) = obj.setlgdNameList(obj.workdirList,obj.model);
        case 'subPlotLayout'
          switch obj.calcType
          case 'Timeseries'
            switch obj.plotType
              case {'Parameter','ParaTimeVar','HyperParameter'}
                obj.(propName) = MyFunc.gridXYSubplot(obj.categoryList); 
              otherwise
                obj.(propName) = [length(obj.categoryList),1]; 
            end
          otherwise
            obj.(propName) = MyFunc.gridXYSubplot(obj.categoryList); 
          end
        case {'plotObj','plotObjY'} 
          if isa(obj,'Title'), return; end
          if ~isempty(obj.lagTermList), obj.setPlotObjLagXX; return; end 
          obj.(propName) = obj.setPlotObj(propName); 
        case 'categoryListY'
          obj.(propName) = obj.categoryList;
        case 'periodListY'
          obj.(propName) = obj.periodList;
        otherwise
          error(strcat("Need to set obj.",propName, "!"))
        end
        disp(strcat("obj.",propName, " is now allocated."))
      else
        disp(strcat("obj.",propName, " is already allocated."))
        switch propName
        case 'lgdNameList'
          if size(obj.lgdNameList,2) ~= size(obj.workdirList,2)
            error('size(obj.workdirList,2) ~= size(obj.lgdNameList,2)!'); 
          end
        case 'periodListY'
          if isempty(obj.plotTypeY) && ~isequal(obj.periodList,obj.periodListY)
            obj.(propName) = obj.periodList;
          end
        end
      end
    end

%% 
    function clearProp(obj,propList)
      for prop = propList
        obj.(prop) = [];
      end
    end

%% 
    function setCalcType(obj,val)
      switch class(val)
      case 'string'
        obj.calcType = val;
      otherwise
        error(strcat("Class of calcType is invalid: ",class(val)))
      end
    end

%%
    function setDataTypeList(obj,val)
      switch class(val)
      case 'string'
        obj.dataTypeList = val;
      otherwise
        error(strcat("Class of dataTypeList is invalid: ",class(val)))
      end
    end

%%
    function setStatsNameList(obj,val)
      switch class(val)
      case 'string'
        obj.statsNameList = val;
        for name = obj.statsNameList
          switch name
            case {'ave','prc68','prc95'}
              %valid: do nothing
            otherwise
              error(strcat("Invalid statsNameList: ",name))
          end
        end
      otherwise
        error(strcat("Class of statsNameList is invalid: ",class(val)))
      end
    end

%% 
    function setWorkdirList(obj,val)
      switch class(val)
        case 'string'
          obj.workdirList = val;
        otherwise
          error(strcat("Class of workdirList is invalid: ",class(val)))
      end
    end

%%
    function y = setPlotObj(obj,propName)
      switch propName
        case 'plotObj'
          plotTypee = "plotType";  categoryListt = "categoryList"; 
        case 'plotObjY'
          plotTypee = "plotTypeY"; categoryListt = "categoryListY";  
      end

      for i = 1:length(obj.workdirList)
        obj.model = obj.model.setup(obj.workdirList(i)); 
        if i == 1
          y = MyCalc.SetPlotObj(obj.(plotTypee),obj.model,obj.(categoryListt));
          if length(obj.workdirList) > 1, y(length(obj.workdirList)) = y; end
        else
          y(i) = MyCalc.SetPlotObj(obj.(plotTypee),obj.model,obj.(categoryListt));
        end
      end
    end

%% 
    function setPlotObjLagXX(obj)
      if length(obj.workdirList) ~= length(obj.lagTermList)
        error('length(obj.workdirList) ~= length(obj.lagTermList)')
      end
      workdirList_org = obj.workdirList;

      for i = 1:length(obj.lagTermList)
        obj.setWorkdirList(workdirList_org(i));
        obj.model.lagFilterTerm = obj.lagTermList(i);
        if i == 1
          obj.plotObj = obj.setPlotObj("plotObj");
          obj.plotObj(length(length(obj.lagTermList))) = obj.plotObj;
        else
          obj.plotObj(i) = obj.setPlotObj("plotObj");
        end
      end
      obj.setWorkdirList(workdirList_org);
    end

%% 
    function checkErrorHistMean(obj,plotObj,nameList)
      if isa(obj,'Gplotmatrix')
        error('You need to use plotMean@Gplotmatrix to mean plotObj');
      end

      if ~strcmp(class(plotObj(1)),class(plotObj(2)))
        error('Different class of plotObj(1) and plotObj(2)"'); 
      end

      if ~isequal(plotObj(1).tspan,plotObj(2).tspan)
        error('Need to set same tspan object'); 
      end
      
      for name = nameList
        if isempty(plotObj(1).(name)) || isempty(plotObj(2).(name)), continue; end
        if size(plotObj(1).(name),1) ~= size(plotObj(2).(name),1) 
          error(strcat('Data length in obj.',name, 'is not consistent!'));
        end
      end
    end
  end

%%
  methods(Abstract)
    setupXYLabelPrep(obj)
    plotDataCategory(obj)
    postDo(obj,varargin)
  end

%%
  methods(Static)
    function lgdNameList = setlgdNameList(workdirList,model)
      lgdNameList = strings(1,length(workdirList));
      for iDir = 1:length(workdirList)
        lgdNameList(iDir) = model.setup(workdirList(iDir)).Info.inventory;
      end
    end

%%
    function y = plotObj3D(plotObj,name)
      y = zeros([size(plotObj(1).(name)),length(plotObj)]); 
      for i = 1:length(plotObj)
        y(:,:,i) = plotObj(i).(name);
      end
    end
  end

end
%
classdef Data < Plot 
  properties
  end

  properties (SetAccess = protected)
    calcType      
    dataTypeList  = ["Prior","Posterior"]
    statsNameList = ["ave"]
  end

  methods
%% 
    function WriteStatics(obj,varargin) 
      obj.addProp(varargin{:});

      objBest = MyCalc.SetPlotObj(obj.plotType,obj.model,obj.categoryList);
      obj.workdir = obj.workdirList(1);

      fname = strcat(obj.workdir,'/outStatics_',erase(obj.plotType,'.'),'_',...
        obj.model.emsCH4File);
      fid = fopen(fname,'w');
      
      fprintf(fid,"%s\n",obj.workdir);
      fprintf(fid,"%s\t%s\n","Para file:",obj.model.inputFile);
      fprintf(fid,"%s\t%s\n","ECH4 file:",obj.model.emsCH4File);
      fprintf(fid,"%s\t%s\n","objName:",obj.plotType);
      fprintf(fid,"%-8s\t%-6s\t%-6s\t%8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\t%-8s\n",...
        "Category","yr_st","yr_end","Ave.","Stdev","Median","Min","Max","2.5Per","97.5Per","5Per","95Per","16Per","84Per","25Per","75Per");
      
      obj.calcWriteStatsLoop(fid,objBest,obj.categoryList,obj.periodList)
      
      fclose(fid);
    end

%%
    function WriteStaticsMultiScenarioHistMean(obj,varargin)
      obj.addProp(varargin{:});
      
      disp(strcat('END WriteStaticsMultiScenarioHistMean ',obj.plotType))      

      inventoryList = reverse(extractBefore(reverse(obj.workdirList),"_")); 
      obj.workdir = obj.workdirList(end);
      fname = strcat(obj.workdir,'/outStatics_',erase(obj.plotType,'.'),'_HistMean_',...
        strjoin(inventoryList,"_"),"_",obj.calcType,'.txt'); 
      fid = fopen(fname,'w');
      
      for i = 1:length(obj.workdirList)
        info = ReadModelParameter(obj.workdirList(i)); 
        fprintf(fid,"%s\n",info.workdir);
        fprintf(fid,"%s\t%s\n","Para file:",info.atmosTargetFile);
        fprintf(fid,"%s\t%s\n","ECH4 file:",info.inventory+".txt");
        fprintf(fid,"%s\t%s\n","objName:",obj.plotType);
      end 

      stats = obj.statsHistMean(obj.plotObj,obj.categoryList,obj.periodList); 

      if strcmp(obj.calcType,'PeriodMean')
        fprintf(fid,"%-8s\t%-6s\t%-6s\t%9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\n",...
        "Category","yr_st","yr_end","Ave.","Median","prc16","prc84","prc2.5","prc97.5","cdf50","cdf16","cdf84","cdf2.5","cdf97.5");
        for iRow = 1:size(obj.periodList{1},1)
        for name = obj.categoryList
          tspanTarget = obj.periodList{1}(iRow,:);
          fprintf(fid,"%-8s\t%6d\t%6d\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\n",...
                    name,tspanTarget,stats.(name).ave(iRow),stats.(name).med(iRow),...
                    stats.(name).prc68(iRow,:),stats.(name).prc95(iRow,:),stats.(name).cdf50(iRow),stats.(name).cdf68(iRow,:),stats.(name).cdf95(iRow,:));
        end
        end
      elseif strcmp(obj.calcType,'PeriodDiff')
        fprintf(fid,"%-8s\t%-6s\t%-6s\t%-6s\t%-6s\t%9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\t%-9s\n",...
        "Category","yr1_st","yr1_end","yr2_st","yr2_end","Ave.","Median","prc16","prc84","prc2.5","prc97.5","cdf50","cdf16","cdf84","cdf2.5","cdf97.5");
        peropdListHorz = horzcat(obj.periodList{:});
        tspanTarget1 = reshape(peropdListHorz(1,:).',2,[]).'; %1st period
        tspanTarget2 = reshape(peropdListHorz(2,:).',2,[]).'; %2nd period
        for iRow = 1:length(obj.periodList)
        for name = obj.categoryList
          tspan1 = tspanTarget1(iRow,:);
          tspan2 = tspanTarget2(iRow,:);
          fprintf(fid,"%-8s\t%6d\t%6d\t%6d\t%6d\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\n",...
                    name,tspan1,tspan2,stats.(name).ave(iRow),stats.(name).med(iRow),...
                    stats.(name).prc68(iRow,:),stats.(name).prc95(iRow,:),stats.(name).cdf50(iRow),stats.(name).cdf68(iRow,:),stats.(name).cdf95(iRow,:));
        end
        end
      end

      fclose(fid);
      disp(strcat('END WriteStaticsMultiScenarioHistMean ',obj.plotType))
    end

%% abstract methods not used for Data.m
    function postDo(obj)
    end
    function plotDataCategory(obj)
    end
    function setupXYLabelPrep(obj)
    end    
  end

  methods(Static)
    function calcWriteStatsLoop(fid,objBest,fieldList,periodList)
      if isa(periodList,'cell')
        periodList = periodList{1};
      end

      for name = fieldList
        for iRow = 1:size(periodList,1)
          tspanTarget = periodList(iRow,1:2);
          objTarget   = MyCalc.PeriodMean(objBest.tspan,objBest.(name),tspanTarget);
          Data.WriteStats(fid,name,tspanTarget,objTarget);          
        end
      end
    end

%% 
    function WriteStats(fid,name,tspan,objTarget)
      ave = mean(objTarget,2);
      sd  = std(objTarget,0,2);
      med = median(objTarget,2);
      MinMax   = prctile(objTarget,[0 100]);
      p95tile  = prctile(objTarget,[2.5 97.5],2);
      p90tile  = prctile(objTarget,[5  95],2);
      p68tile  = prctile(objTarget,[16 84],2);
      p50tile  = prctile(objTarget,[25 75],2);
      fprintf(fid,"%-8s\t%6d\t%6d\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\t%9.5f\n",...
              name,tspan,ave,sd,med,MinMax,p95tile,p90tile,p68tile,p50tile); 
    end

%%
    function plotObj = ReadOutStaticsMulti(workdir,plotType,calcType,periodList)
      if length(periodList) ~= 1, error("length(periodList) ~= 1"); end 

      info = ReadModelParameter(workdir);
      inventoryList = reverse(extractBefore(reverse(info.workdirList),"_"));

      fname = strcat(workdir,'/outStatics_',erase(plotType,'.'),'_HistMean_', ...
        strjoin(inventoryList,"_"),"_",calcType,'.txt');
      [fid,errmsg] = fopen(fname,'r'); 

      if ~isempty(errmsg) && strcmp("PeriodMean",calcType)
        [fid,errmsg] = fopen(erase(fname,strcat("_",calcType)),'r'); 
      end

      if ~isempty(errmsg), disp(fname); error(errmsg); end

      textscan(fid,'%s',4*length(info.workdirList)+1,'delimiter','\n');

      plotObj = [];
      while ~feof(fid) 
        name = textscanString(fid,'%s',1,'delimiter','\t');
        switch calcType
          case "PeriodDiff"
            period = cell2mat(textscan(fid,'%f',4,'delimiter','\t'));
          otherwise
            period = cell2mat(textscan(fid,'%f',2,'delimiter','\t'));
        end
    
        if isequal(period,reshape(periodList{:}.',[],1))
          plotObj.(name).ave    = cell2mat(textscan(fid,'%f',1,'delimiter','\t'));
          plotObj.(name).med    = cell2mat(textscan(fid,'%f',1,'delimiter','\t'));
          plotObj.(name).p68    = cell2mat(textscan(fid,'%f',2,'delimiter','\t'));
          plotObj.(name).p95    = cell2mat(textscan(fid,'%f',2,'delimiter','\t'));
          textscan(fid,'%s',1,'delimiter','\n');
        else
          textscan(fid,'%s',1,'delimiter','\n');
        end
      end  
      if isempty(plotObj), error("No corresponding periodList exists!"); end
    end  

  end
end

%%
function y =  textscanString(fid,format,num,~,delimiter)
  y_cell = textscan(fid,format,num,'delimiter',delimiter);
  y      = erase(convertCharsToStrings(y_cell{1})," ");      
end
    
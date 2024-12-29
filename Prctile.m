% A subclass of Plot.m. It creates the timeseries of each 
% output (e.g. atmosphere, parameters) with prescribed percentiles. 
% Ryo Fujita 2024

classdef Prctile < Plot 
  properties (SetAccess = protected) 
    calcType       = "Timeseries"
    dataTypeList   
    statsNameList  = ["prc68","prc95","ave"];
  end

  methods
%%
    function obj = Prctile(varargin)
      obj = obj@Plot(varargin{:}); 
      switch obj.plotType
      case {'Parameter','d13CdDKIECKIEDtot','ParaTimeVar','HyperParameter'}
        obj.dataTypeList = "Posterior";
      otherwise
        obj.dataTypeList = ["Prior","Posterior"];
      end
    end

%% 
    function h = postDo(obj) 
      obj.modStyle.plotLegend(obj.subPlotLayout);

      if obj.subPlotLayout(1) == 1 || obj.subPlotLayout(2) > 1
        obj.modStyle.isUseSameAxis = "off"; 
      end

      obj.modStyle.XYLabelTimePost(obj.plotType);
      obj.modStyle.setXlimTimeSubplot(obj.periodList{1},obj.subPlotLayout(2)); 
      obj.modStyle.useSameAxis(obj.categoryList); 
      obj.modStyle.setSizeFig(obj.subPlotLayout); 
      obj.modStyle.setFigTransparency; 
      obj.modStyle.isUseSameAxis = "on";
    end

%%
    function setupXYLabelPrep(obj)
      if ~strcmp(obj.workdir,obj.workdirList(1)), return; end
      grid on; hold on
      xlabel('Year'); 
      switch obj.plotType
        case 'Atmosphere'
          ylabel(MyFunc.nameLabel(obj.currentCategory,obj.plotType)); 
          if strcmp(obj.currentCategory,obj.categoryList(1))
            MyFunc.setPlotAtmosTarget(obj.categoryList,obj.subPlotLayout,obj.model,obj.periodList{1});
          end 
        case {'ECH4','SourceFraction','ECH4life','ECH4NPP'} 
          MyFunc.setTextCategory(obj.currentCategory); 
        otherwise
          ylabel(MyFunc.nameLabel(obj.currentCategory,obj.plotType)); 
      end
    end

%%
    function obj = plotDataCategory(obj) 
      if isempty(obj.plotObj(obj.iDir).(obj.currentCategory)), return; end
      tscale = TspanScale(obj.plotObj(obj.iDir).tspan,obj.periodList{1});

      for dataType = obj.dataTypeList
        h = obj.plotDataType(tscale.tspan,dataType); 
        iCategory = find(ismember(obj.categoryList,obj.currentCategory));
        obj.modStyle = obj.modStyle.setLegend(h,dataType,obj.inventory,iCategory);
      end

      obj.modStyle.setYAxisLimitCategoryTime(...
        obj.plotType,obj.currentCategory,obj.periodList{1},obj.model);
    end

%%
    function h = plotDataType(obj,tspan,dataType)
      switch dataType
        case "Prior"
          [tspanPrior,dataPrior] = obj.YdataPrior;
          if isempty(dataPrior), h = []; return; end

          tscale = TspanScale(tspanPrior,obj.periodList{1});
          h = plot(tscale.tspan,dataPrior);
          set(h,'Color',obj.modStyle.color(obj.inventory,obj.iDir))
          set(h,'LineWidth',2)
          set(h,'LineStyle',':')
        case "Posterior"
          for statsName = obj.statsNameList
            switch statsName
            case "ave"
              h = plot(tspan,obj.YdataBestCombination(statsName),...
                '-','DisplayName',obj.inventory);
              h.Color = obj.modStyle.color(obj.inventory,obj.iDir);

              if length(obj.dataTypeList) > 1
                h.LineWidth = 2; 
              else
                h.LineWidth = 1.75; 
              end
            case {"prc68","prc95","prc90"}
              h = MyFunc.plotArea(tspan,obj.YdataBestCombination(statsName)); 
              h = obj.modStyle.colorPrcTime(obj.inventory,obj.iDir,statsName,obj.statsNameList,h);
            end
          end
      end
    end

%%
    function [tspan,y] = YdataPrior(obj)
      switch obj.plotType
        case {"ECH4","ECH4life"} 
          tspan = obj.plotObj(obj.iDir).tspan;
          E_CH4_org = Emission().ReadFileCH4(obj.model,tspan);
          E_CH4_org = convertEmsFromScalingFactor(E_CH4_org,obj.model); 
          if sum(ismember(fieldnames(E_CH4_org),obj.currentCategory)) == 0, y = [];
          else, y = E_CH4_org.(obj.currentCategory);
          end
        case "SourceFraction"
          tspan = obj.plotObj(obj.iDir).tspan;
          E_CH4_org = Emission().ReadFileCH4(obj.model,tspan);
          E_CH4_org = convertEmsFromScalingFactor(E_CH4_org,obj.model);  
          y = E_CH4_org.(obj.currentCategory)./E_CH4_org.tot .* 100;
        case "Atmosphere" 
          obj.model.randType = 'Prior'; %obj.model.status   = 10; 
          obj.model.numCase = 1; 
          obj.model = obj.model.mainWithoutFiltering; 
          obj.model.atmos.CH4 = obj.model.atmos.CH4./MassConstants.factor; 
          y = obj.model.atmos.(obj.currentCategory);
          tspan = obj.model.atmos.tspan;
        otherwise
          error(strcat("Prior data in ",obj.plotType," is not defined!"))
      end
    end

%%
    function y = YdataBestCombination(obj,statsName)      
      if ~isa(obj.plotObj(obj.iDir).(obj.currentCategory),'double') && ...
         any(contains(statsName,string(fieldnames(obj.plotObj(obj.iDir).(obj.currentCategory)))))  
        y = obj.plotObj(obj.iDir).(obj.currentCategory).(statsName);
      else
        switch statsName
          case "ave"
            y = mean(obj.plotObj(obj.iDir).(obj.currentCategory),2);
          case {"prc68"}
            y = prctile(obj.plotObj(obj.iDir).(obj.currentCategory),[16 84],2);
          case {"prc95"}
            y = prctile(obj.plotObj(obj.iDir).(obj.currentCategory),[2.5 97.5],2);
          case {"prc90"} 
            y = prctile(obj.plotObj(obj.iDir).(obj.currentCategory),[5 95],2);
        end
      end
    end

%%
    function setStatsNameList(obj,statsNameList)
      if isa(statsNameList,"string")
        obj.statsNameList = statsNameList;
      end
    end

  end
end

%%
function E_CH4_org = convertEmsFromScalingFactor(E_CH4_org,model)
  if contains("fbb",model.nameList.para) 
    fbb = mean(vertcat(vertcat(model.paraTargetList.fbb).def));
    Ebb_org = E_CH4_org.bb;
    E_CH4_org.bb = Ebb_org.*fbb;
    E_CH4_org.anth = E_CH4_org.anth - Ebb_org + E_CH4_org.bb;
    E_CH4_org.tot = E_CH4_org.tot - Ebb_org + E_CH4_org.bb;
  end
end
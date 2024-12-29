% A subclass of Plot.m. It creates a histogram of each 
% output (e.g. atmosphere, parameters). 
% Ryo Fujita 2024

classdef Histogram < Plot
  properties(Constant)
    NumBins = 20; %default
  end

  properties (SetAccess = protected)
    calcType       = "PeriodMean";
    dataTypeList   = "Posterior";
    statsNameList  = "ave";
  end

  methods
    function obj = plotDataCategory(obj)      
      if size(obj.dataTypeList,2) > 1
        error('Currently size(obj.dataTypeList,2) in Histogram needs to be one!'); 
      end

      set(gca,'FontSize',12); box on

      y = obj.plotObj(obj.iDir).(obj.currentCategory);
      
      if isempty(y)
        disp(strcat("obj.plotObj(obj.iDir).",obj.currentCategory," is empty!")); 
        return; 
      end 

      if isa(y,'struct')
        for i = 1:size(y.BinEdges,1)
          h = histogram('BinEdges',y.BinEdges(i,:),'BinCounts',y.Values(i,:), ...
            'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
        end
      else   
        y = MyCalc.convertCalcType(obj.calcType,obj.plotObj(obj.iDir), ...
          obj.currentCategory,obj.periodList).(obj.currentCategory);
        for i = 1:size(y,1)
          h = histogram(y(i,:),'Normalization','probability','DisplayStyle', ...
            'stairs','LineWidth',2); hold on
          h.BinEdges = MyFunc.EdgePara(class(obj.plotObj(obj.iDir)), ...
            obj.model,obj.currentCategory);
          if size(obj.periodList{1},1) > 1 && length(obj.workdirList) == 1 ...
            && size(obj.lagTermList,1) <= 1 
            c = MyFunc.hsv_v2(size(y,1));
            h.EdgeColor = c(i,:);  
          else 
            h.EdgeColor = obj.modStyle.color(obj.inventory,obj.iDir);
          end
        end
      end      

      iCategory    = find(ismember(obj.categoryList,obj.currentCategory));
      obj.modStyle = obj.modStyle.setLegend(h,obj.dataTypeList,obj.inventory,iCategory);
    end

%% 
    function setupXYLabelPrep(obj)
      if ~strcmp(obj.workdir,obj.workdirList(1)), return; end
      title(MyFunc.nameLabel(obj.currentCategory,obj.plotType)); hold on
    end

%% 
    function postDo(obj)
      obj.plotLegend;
      obj.modStyle.setSizeFig(obj.subPlotLayout); 
    end

%%
    function updateHistCommonLim(obj)
      if length(obj.plotObj) == 1 || ~isa(obj.plotObj(1).(obj.categoryList(1)),'struct'), return; end 

      for iCat = 1:length(obj.categoryList)
        name = obj.categoryList(iCat);
        numBinsRow = zeros(1,length(obj.plotObj)); 
        binLimitsRow = zeros(length(obj.plotObj),2); 

        for i = 1:length(obj.plotObj)
          numBinsRow(i)     = obj.plotObj(i).(name).NumBins; 
          binLimitsRow(i,:) = obj.plotObj(i).(name).BinLimits; 
        end

        if (unique(numBinsRow) == obj.NumBins ...
        && (length(unique(binLimitsRow(:,1))) == 1 ...
        && length(unique(binLimitsRow(:,2))) == 1)) ...
        || ~any(binLimitsRow,'all')
          continue
        end

        [binLimMin, binLimMax] = bounds(binLimitsRow,'all'); 
        binEdges  = binLimMin:(binLimMax-binLimMin)/obj.NumBins:binLimMax;

        for i = 1:length(obj.plotObj)
          values = zeros(1,obj.NumBins);
          for j = 1:obj.NumBins
            for k = 2:obj.plotObj(i).(name).NumBins+1
              if obj.plotObj(i).(name).BinEdges(k) > binEdges(j) ...
              && obj.plotObj(i).(name).BinEdges(k) <= binEdges(j+1) 
                values(j) = values(j) + obj.plotObj(i).(name).Values(k-1);
              end
            end
          end

          if round(sum(values,'all'),4) ~= 1
            error('Unproper solution in HistogramMean'); 
          else
            obj.plotObj(i).(name).NumBins   = obj.NumBins;
            obj.plotObj(i).(name).BinLimits = [binLimMin, binLimMax];
            obj.plotObj(i).(name).BinEdges  = binEdges;
            obj.plotObj(i).(name).Values    = values;
          end
        end
      end
    end

%%
    function plotLegend(obj)
      if strcmp(obj.calcType,'PeriodMean') && size(obj.periodList{1},1) > 1 ...
      && size(obj.workdirList,1) == 1 && size(obj.lagTermList,1) <= 1 
        for i = 1:size(obj.periodList{1},1)
          obj.modStyle.lgd.String(i) = {strcat(num2str(obj.periodList{1}(i,1)), ...
          '-',num2str(obj.periodList{1}(i,2)))};
        end

        obj.modStyle.lgd.Location    = 'northwest';  
        obj.modStyle.lgd.Orientation = 'horizontal'; 
        obj.modStyle.lgd.FontSize    = 8;
        obj.modStyle.lgd.Box = 'off'; 

        if obj.subPlotLayout(1) == 3 && obj.subPlotLayout(2) == 4
          obj.modStyle.lgd.NumColumns = 4;
          obj.modStyle.lgd.Position = [0.52 0.13 0.9062 0.1321];
        else
          obj.modStyle.lgd.Position = [0.3    0.0056    0.3045    0.0333]; 
          if length(obj.modStyle.lgd.String) > 5 
            if obj.subPlotLayout(1) == 1 && obj.subPlotLayout(2) == 1
              obj.modStyle.lgd.NumColumns  = 2; 
              obj.modStyle.lgd.Location = 'best';
            else
              obj.modStyle.lgd.NumColumns  = 5; 
            end
          end 
        end
      else
        obj.modStyle.plotLegend(obj.subPlotLayout);
      end
    end
  
  end
end
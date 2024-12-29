%

classdef BarPlot < Plot
  properties
    plotStyle = "errorbar" %'bar', 'boxplot', 'errorbar' etc... 
    xval = 0.5;
    plotPriorStyle = 'yes'; %'no', 'yes', 'single'
    xTickLabel 
    setXLabel = 'yes'; %'yes' or 'no'
  end

%%
  properties (SetAccess = protected)
    calcType       = "PeriodMean"
    dataTypeList   = ["Prior","Posterior"]
    statsNameList  = ["ave","prc68"];
  end

%%
  methods 
    function setStatsNameList(obj,val)
      if length(val) < 2
        error('statsNameList length in BarPlot needs to be more than 2'); 
      end
  
      setStatsNameList@Plot(obj,val);
    end

%%
    function obj = plotDataCategory(obj)
      if strcmp(obj.plotType,'Parameter') || strcmp(obj.flagScenarioMean,'y')
        obj.plotPriorStyle = 'single';
      end

      for dataType = obj.dataTypeList
        obj.currentDataType = dataType;
        if isempty(obj.plotObj(obj.iDir).(obj.currentCategory))
          continue; 
        end

        switch obj.currentDataType
          case "Prior"
            if strcmp(obj.plotPriorStyle,'no') || (strcmp(obj.plotPriorStyle,'single') ...
                && ~strcmp(obj.workdir,obj.workdirList(1)))
              continue;
            end

            y = obj.YdataPrior(obj.plotType);
            if size(y,2) == 1
              for i = 1:size(y,1)
                h = plot(y(i),'_','MarkerEdgeColor','k','MarkerSize',8,'LineWidth',0.8);
                modifyYLimLinePlot(y(i));
              end

            elseif size(y,2) == 3
              for i = 1:size(y,1)
                h = errorbar(1,y(i,1),y(i,2),y(i,3),'k_','MarkerSize',8,'CapSize',8, ...
                  'MarkerEdgeColor','k','MarkerFaceColor','k');
                hold on
                modifyYLimErrorBar(y(i,:));
              end
            end

          case "Posterior"
            switch obj.plotStyle
              case 'bar'
                h = bar(obj.YdataBestBar);
                hold on
              case 'errorbar'
                y = obj.YdataBestErrorBar;
                for i = 1:size(y,1)
                  h = errorbar(1,y(i,1),y(i,2),y(i,3),'ks','MarkerSize',8,'CapSize',8, ...
                    'MarkerEdgeColor',obj.modStyle.color(obj.inventory,obj.iDir), ...
                    'MarkerFaceColor',obj.modStyle.color(obj.inventory,obj.iDir));
                  hold on

                  if ismember("prc68",obj.statsNameList)
                    h.Color = obj.modStyle.color(obj.inventory,obj.iDir);
                    h.LineWidth = 2;
                    h.CapSize = 0;
                  end

                  modifyYLimErrorBar(y(i,:));
                end
            end
        end

        set(gca,'FontSize',12);
        iCategory    = find(ismember(obj.categoryList,obj.currentCategory));
        obj.modStyle = obj.modStyle.setLegend(h,dataType,obj.inventory,iCategory);
      end
    end

%%
    function y = YdataPrior(obj,plotType)
      switch plotType
        case "Parameter" 
          if strcmp(obj.plotStyle,'errorbar')       
            y = obj.getPriorErrorBar("Parameter",obj.currentCategory);
          end
          return

        case "ECH4"
          if strcmp(obj.flagScenarioMean,'y')
            for i = 1:length(obj.model.Info.workdirList)
              obj.model = obj.model.setup(obj.model.Info.workdirList(i));
              E_CH4_org = Emission().ReadFileCH4(obj.model,obj.plotObj(obj.iDir).tspan);
              y.tspan   = E_CH4_org.tspan; 
              
              if i == 1, y_tmp = zeros(length(y.tspan),2); end

              y_tmp(:,i)= E_CH4_org.(obj.currentCategory);
            end
            y.(obj.currentCategory) = mean(y_tmp,2);
          else
            E_CH4_org = Emission().ReadFileCH4(obj.model,obj.plotObj(obj.iDir).tspan);
            y.tspan   = E_CH4_org.tspan;
            y.(obj.currentCategory) = E_CH4_org.(obj.currentCategory);
          end

        case "SourceFraction"
          if strcmp(obj.flagScenarioMean,'y')
            for i = 1:length(obj.model.Info.workdirList)
              obj.model = obj.model.setup(obj.model.Info.workdirList(i));
              E_CH4_org = Emission().ReadFileCH4(obj.model,obj.plotObj(obj.iDir).tspan);
              y.tspan   = E_CH4_org.tspan; 

              if i == 1, y_tmp = zeros(length(y.tspan),2); end
              
              y_tmp(:,i)= E_CH4_org.(obj.currentCategory)./E_CH4_org.tot;
            end
            y.(obj.currentCategory) = mean(y_tmp,2);
          else
            E_CH4_org = Emission().ReadFileCH4(obj.model,obj.plotObj(obj.iDir).tspan);
            y.tspan   = E_CH4_org.tspan;
            y.(obj.currentCategory) = E_CH4_org.(obj.currentCategory)./E_CH4_org.tot .* 100;
          end

        case "Atmosphere" %2021.2.26
          obj.model.randType = 'Prior'; %2021.5.10 TBD
          obj.model.status   = 10; %2021.2.26 TBD
          if strcmp(obj.flagScenarioMean,'y')
            for i = 1:length(obj.model.Info.workdirList)
              y.tspan    = Atmosphere(obj.model).tspan; 

              if i == 1, y_tmp = zeros(length(y.tspan),2); end
              
              y_tmp(:,i) = Atmosphere(obj.model).(obj.currentCategory); 
            end
            y.(obj.currentCategory) = mean(y_tmp,2);
          else
            y.tspan = Atmosphere(obj.model).tspan;
            y.(obj.currentCategory) = Atmosphere(obj.model).(obj.currentCategory); 
          end

        otherwise
          if sum(ismember(obj.model.nameList.para,obj.currentCategory)) == 1 
            y = obj.YdataPrior("Parameter");
          elseif sum(ismember(string(fieldnames(SourceWorking)).',obj.currentCategory)) == 1 
            y = obj.YdataPrior("ECH4");
          elseif strcmp(obj.currentCategory,"life") 
            switch obj.calcType 
              case "PeriodMean", y = [9.1,9.1-7.7,9.5-9.1]; %floss +/-10%
              case "PeriodDiff", y = [0 0 0];
              otherwise, error('Not yet prepoared!')
            end
          else
            error('No such a property name both in SourceWorking or ParameterSourceSink!');
          end
          return
      end

      y = MyCalc.convertCalcType(obj.calcType,y,obj.currentCategory,obj.periodList{1});

      if strcmp(obj.calcType,'PeriodMean') && strcmp(obj.plotStyle,'errorbar')
        y = obj.getPriorErrorBar(y.(obj.currentCategory));
      else
        y = y.(obj.currentCategory);
      end
    end

%%
    function y = getPriorErrorBar(obj,varargin)
      if ~isempty(varargin) && nargin == 2
        y_org        = varargin{1};
        plotTypeName = obj.plotType;
        catName      = obj.currentCategory;
      elseif ~isempty(varargin) && nargin == 3 
        plotTypeName = varargin{1};
        catName      = varargin{2};
      else
        plotTypeName = obj.plotType;
        catName      = obj.currentCategory;
      end

      if strcmp(catName,"life") && ~ismember("life",obj.model.nameList.para)
        switch obj.calcType
          case "PeriodMean", y = [9.1,9.1-7.7,9.5-9.1]; %floss +/-10%
          case "PeriodDiff", y = [0 0 0];
          otherwise, error('Not yet prepoared!')
        end
        return
      end

      switch plotTypeName
        case 'Parameter'
          if strcmp(obj.calcType,'PeriodMean')
            y = [obj.model.paraTargetList(1).(catName).def, ...
              abs(obj.model.paraTargetList(1).(catName).def ...
              - obj.model.paraTargetList(1).(catName).min),...
              abs(obj.model.paraTargetList(1).(catName).def ...
              - obj.model.paraTargetList(1).(catName).max)]; 
          else 
            y = 0;
          end
        case {'ECH4','ECH4life'}
          switch catName
            case {'anth_bio','natr_bio','anth_ff','bb'}
              para_error = obj.getPriorErrorBar("Parameter",strcat('f',catName));
              y = y_org .* para_error;
            case 'geo'
              y = obj.getPriorErrorBar("Parameter",strcat('E',catName));
            otherwise
              disp(strcat(plotTypeName,".",catName,": prior errorbar is not yet prepared!"))
              y = y_org;
          end
        case 'SourceFraction'
          disp(strcat(plotTypeName,".",catName,": prior errorbar is not yet prepared!"))
          y = y_org;
      end
    end

%%
    function y = YdataBestBar(obj)
      if isa(obj.plotObj(obj.iDir).(obj.currentCategory),'struct') 
        pObj = obj.plotObj(obj.iDir).(obj.currentCategory);
        y = pObj.(obj.currentCategory).ave; 
      else
        pObj = MyCalc.convertCalcType(obj.calcType,obj.plotObj(obj.iDir), ...
          obj.currentCategory,obj.periodList);
        y = mean(pObj.(obj.currentCategory),2); 
      end
    end

%%
    function y = YdataBestErrorBar(obj)    
      if isa(obj.plotObj(obj.iDir).(obj.currentCategory),'struct')
        pObj = obj.plotObj(obj.iDir).(obj.currentCategory);
        y    = pObj.(obj.statsNameList(1)); 
        y(2) = abs(pObj.ave-pObj.(obj.statsNameList(2))(1));
        y(3) = abs(pObj.ave-pObj.(obj.statsNameList(2))(2));

      else
        pObj = MyCalc.convertCalcType(obj.calcType,obj.plotObj(obj.iDir), ...
          obj.currentCategory,obj.periodList);

        switch obj.statsNameList(1)
        case 'ave'
          y = mean(pObj.(obj.currentCategory),2);
        case 'med'
          y = median(pObj.(obj.currentCategory),2);
        end

        prcXX = MyCalc.(obj.statsNameList(2))(pObj.(obj.currentCategory));
        y(:,2) = abs(y(:,1)-prcXX(:,1));
        y(:,3) = abs(y(:,1)-prcXX(:,2));
      end
    end

%%
    function postDo(obj)
      obj.modStyle.plotLegend(obj.subPlotLayout);
      obj.modifyPlotStyle; %only BarPlot
      obj.modStyle.setSizeFig(obj.subPlotLayout);
    end

%%
    function modifyPlotStyle(obj)
      fig = gcf;
      ax  = findobj(fig,'Type','Axes');

      for i = 1:length(ax)
        fig.CurrentAxes = ax(i);
        box on;
        obj.currentCategory = obj.categoryList(length(ax)-i+1);

        pObj = ax(i).Children;
        for j = 1:length(pObj)     
          pObj(j).XData = length(pObj) - j + 1;  
        end
        
        ax(i).XTick = 1:1:length(pObj);

        if size(obj.periodList{1},1) > 1 && length(obj.workdirList) == 1
          if strcmp(obj.plotPriorStyle,'single') ...
            && length(pObj) == 2*size(obj.periodList{1},1)
            obj.plotPriorStyle = 'singleMultiPeriod'; 
          end

          if strcmp(obj.calcType,'PeriodDiff') 
            lgdList = strings(1,size(obj.periodList,1));
            for iprd = 1:size(obj.periodList,1)
              lgdList(iprd) = ...
              [strcat(extractAfter(num2str(obj.periodList{1}(2,1)),2), ...
              '-',extractAfter(num2str(obj.periodList{1}(2,2)),2)),...
              ' vs ',strcat(extractAfter(num2str(obj.periodList{1}(1,1)),2), ...
              '-',extractAfter(num2str(obj.periodList{1}(1,2)),2))];
            end
          else
            lgdList = strings(1,size(obj.periodList{1},1));
            for iprd = 1:size(obj.periodList{1},1)
              lgdList(iprd) = ...
              {strcat(num2str(obj.periodList{1}(iprd,1)), ...
              '-',num2str(obj.periodList{1}(iprd,2)))};
            end
          end
        else
          lgdList = obj.lgdNameList;
        end

        switch obj.plotPriorStyle
          case 'single'
            ax(i).XTickLabel = ['Prior' lgdList];
          case 'no'
            ax(i).XTickLabel = lgdList;
          case {'yes','singleMultiPeriod'}
            ax(i).XTickLabel = {};
            if strcmp(obj.setXLabel,'yes')
              for lgd = lgdList
                ax(i).XTickLabel = [ax(i).XTickLabel; strcat(["Prior-","Post-"],lgd).'];
              end
            else
              for lgd = lgdList 
                ax(i).XTickLabel = [ax(i).XTickLabel; ["Prior","Post"].'];
              end
            end
        end
        
        if obj.subPlotLayout(2) >= 4 || size(char(ax(i).XTickLabel),2) > 8
          ax(i).XTickLabelRotation = 45; 
        end 
        
        if size(char(ax(i).XTickLabel),2) > 10
          ax(i).FontSize = 10;
        elseif size(char(ax(i).XTickLabel),2) > 8
          ax(i).FontSize = 11;
        else
          ax(i).FontSize = 12;
        end

        if length(pObj) <= 5
          ax(i).XLim(1) = ax(i).XLim(1)-0.5;
          ax(i).XLim(2) = ax(i).XLim(2)+0.5;
        else
          ax(i).XLim(2) = ax(i).XLim(2)+1;
        end
      end
    end

%% 
    function setupXYLabelPrep(obj)
      if ~strcmp(obj.workdir,obj.workdirList(1)), return; end
      ylabel(MyFunc.nameLabel(obj.currentCategory,obj.plotType)); hold on
    end
  end
end

%%
function modifyYLimErrorBar(y)
  ax = gca; 
  ylim_org = ax.YLim;  
  
  if y(1)-y(2) <= ylim_org(1)
    ax.YLim(1) = y(1)-y(2); 
  end
  
  if y(1)+y(3) >= ylim_org(2)
    ax.YLim(2) = y(1)+y(3); 
  end

  ylimRange = ax.YLim(2) - ax.YLim(1);

  if y(1)-y(2) <= ylim_org(1)
    ax.YLim(1) = ax.YLim(1) - ylimRange/10; 
  end

  if y(1)+y(3) >= ylim_org(2)
    ax.YLim(2) = ax.YLim(2) + ylimRange/10; 
  end
end

%%
function modifyYLimLinePlot(y) 
  ax = gca;
  ax.YLim(1) = ax.YLim(1) - y*0.05; 
  ax.YLim(2) = ax.YLim(2) + y*0.05;
end
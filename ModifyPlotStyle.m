% Controls the figure sizes, axis styles, colors, and legends.
% Ryo Fujita 2024
classdef ModifyPlotStyle
  properties 
    isUseSameAxis     = "on";
    setYAxisLocation  = 'alt2'; %'alt2' or 'left'
    setSeparateYAxis  = 'no'; %'yes' or 'no'
    setTransparency   = 'no'; %'yes' or 'no'
    lgd
    model
  end

  methods
    function obj = ModifyPlotStyle(varargin)
      addpath samexaxis
      if nargin > 0
        for k = 1:nargin      
          switch k 
            case 1, obj.model = varargin{1};
          end
        end
      end
    end

%% 
    function obj = setLegend(obj,h,dataType,dataPlotName,iCategory)
      if isempty(h), return; end 
      if isempty(obj.lgd)
        obj.lgd = legend; 
        obj.lgd.UserData.handle = []; 
        obj.lgd.UserData.name = [];
      end

      if iCategory == 1
        dataPlotName = char(dataPlotName);
        if strcmp(dataPlotName(1:1),'_')
          dataPlotName = dataPlotName(2:length(dataPlotName)); 
        end

        switch dataType
          case "Prior"
            hObj = h(1);  
            suffix = strcat(string(dataPlotName),"-Prior");
          case "Posterior"
            if size(h,2) == 1
              hObj = h(1);
            else
              hObj = h(2);
            end
            suffix = strcat(string(dataPlotName),"-Post");
        end
        
        obj.lgd.UserData.handle  = [obj.lgd.UserData.handle,hObj];
        hName = suffix;
        obj.lgd.UserData.name    = [obj.lgd.UserData.name,hName];
      end
    end
 
%% 
    function obj = plotLegend(obj,subPlotLayout)
      if isempty(obj.lgd)
        return
      else
        obj.lgd = legend(obj.lgd.UserData.handle,obj.lgd.UserData.name);
      end

      if subPlotLayout(1) > 1 && subPlotLayout(2) == 1
        obj.lgd.NumColumns  = 2;
        %obj.lgd.Orientation = 'vertical';
        obj.lgd.Orientation = 'horizontal';
        obj.lgd.Location    = 'southeast';
      else
        if length(obj.lgd.String) > 5, obj.lgd.NumColumns  = 5; end 
        obj.lgd.Location    = 'northwest';  
        obj.lgd.Orientation = 'horizontal'; 
        obj.lgd.Position = [0.3    0.0056    0.3045    0.0333]; 
      end

      obj.lgd.FontSize = 8;
      obj.lgd.Box = 'off'; 
    end

%% 
    function XYLabelTimePost(obj,plotType)
      fig = gcf;
      if strcmp(obj.isUseSameAxis,"on") && ~contains(class(fig.Children),"TiledChartLayout")
        samexaxis([],'YAxisLocation',obj.setYAxisLocation,'XMinorTick','on', ...
          'Join','YTickAntiClash','YLabelDistance',0.5);
      end

      switch plotType
        case 'ECH4'
          an = annotation('textarrow','Position',[0.03 0.6 0 0], ...
            'String','CH_{4} Emission (Tg/yr)', 'FontSize',13, ...
            'VerticalAlignment','middle','HeadStyle','none','TextRotation',90);
            pause(0.1);
        case 'SourceFraction'
          an = annotation('textarrow','Position',[0.03 0.6 0 0], ...
            'String','Source Fraction (%)', 'FontSize',13, ...
            'VerticalAlignment','middle','HeadStyle','none','TextRotation',90);
            pause(0.1);
      end
    end

%% 
    function h = colorPrcTime(obj,plotName,iDir,prcName,prctileList,h)
      if size(h,2) == 2
        set(h(2),'FaceColor',obj.color(plotName,iDir),'EdgeColor',obj.color(plotName,iDir)); 
        if strcmp(prcName,"prc95")
          set(h(2),'EdgeAlpha',0.2,'LineStyle','-','LineWidth',0.5,'FaceAlpha',0.20);
        elseif strcmp(prcName,"prc68") 
          set(h(2),'EdgeAlpha',0.2,'LineStyle','-','LineWidth',0.5,'FaceAlpha',0.20); 
        elseif strcmp(prcName,"prc90")
          set(h(2),'EdgeAlpha',0.2,'LineStyle','-','LineWidth',0.5,'FaceAlpha',0.20);
        end
      else
        h.Color = obj.color(plotName,iDir);
      end
    end

%% 
    function useSameAxis(obj,categoryList) 
      if ~strcmp(obj.isUseSameAxis,"on"), return; end
      obj.modifyFigureShape(categoryList);
      obj.separateYAxis;
      obj.shrinkXSize;
    end

%% 
    function modifyFigureShape(obj,categoryList)
      fig = gcf; pause(0.3);
      ax = findobj(fig,'Type','Axes');
      for iCategory = 1:length(categoryList)
        fig.CurrentAxes = ax(length(categoryList)-iCategory+1);
        fig.CurrentAxes.Position(3) = fig.CurrentAxes.Position(3)*0.93;
        fig.CurrentAxes.Position(1) = fig.CurrentAxes.Position(1) + fig.CurrentAxes.Position(3)*0.03;
      end
    end

%%
    function separateYAxis(obj)
      if ~strcmp(obj.setSeparateYAxis,"yes"), return; end 
      ax = findobj(gcf,'Type','Axes');
      numRow  = length(ax);
      yoffset = 1 - (ax(numRow).Position(2) + ax(numRow).Position(4));
      for i = numRow:-1:1
        ax(i).Position(2) = ax(i).Position(2) + (yoffset - 0.02)/(numRow-1)*(i-1);
        if i ~= 1, ax(i).XLabel.String = ''; end
      end
    end 

%% 
    function shrinkXSize(obj)
      if ~strcmp(obj.setYAxisLocation,"yes"), return; end 
      ax = findobj(gcf,'Type','Axes');
      fig = gcf; pause(0.3);
      numRow  = length(ax);
      for i = 1:numRow
        fig.CurrentAxes = ax(i);
        fig.CurrentAxes.Position(3) = fig.CurrentAxes.Position(3)*0.97;
        if mod(numRow-i+1,2) == 0, fig.CurrentAxes.YLabel.Position(1) = 1.05; end
      end
    end

%%
    function setSizeFig(obj,subPlotLayout)
      fig = gcf; pause(0.3);
  
      % convert from inch to cm for ppt layout
      fig.PaperUnits = 'centimeters';

      if subPlotLayout(1) == 1 && subPlotLayout(2) == 1
        fig.CurrentAxes.Position(4) = fig.CurrentAxes.Position(4)*0.95; 
        fig.CurrentAxes.Position(2) = 0.12;
        fig.PaperPosition = [0 0 26 23];
        fig.PaperOrientation = 'landscape';
      elseif subPlotLayout(2) == 1 && strcmp(obj.isUseSameAxis,"on")
        fig.PaperPosition(1) = 0; fig.PaperPosition(2) = 0; 
        fig.PaperPosition(3) = 15; 
        fig.PaperPosition(4) = 26 +(subPlotLayout(1)-4) * 4.5; 
      elseif subPlotLayout(2) == 1
        fig.PaperPosition(1) = 0; fig.PaperPosition(2) = 0; 
        fig.PaperPosition(3) = 14;
        fig.PaperPosition(4) = 25 +(subPlotLayout(1)-4) * 3;
      elseif subPlotLayout(1) == 1
        fig.PaperPosition = [0 0 33 9]; 
      else
        fig.PaperPosition = [0 0 33 17.3]; %default 
      end

      fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)]; 
    end


%% 
    function setFigTransparency(obj)
      if ~strcmp(obj.setTransparency,'yes'), return; end

      fig = gcf;
      ax = findobj(fig,'Type','Axes');
      
      fig.Color = 'None';

      numRow  = length(ax);
      for i = 1:numRow
        ax(i).Color = 'None';
      end
      fig.InvertHardcopy = 'off';
    end

%%
    function setXlimTimeSubplot(obj,xLim,numCol)
      fig = gcf;
      ax  = findobj(fig.Children,'Type','Axes');

      if obj.model.yrInitial >= 2015 && xLim(1) < 1950, xLim(1) = 2010; end
      if obj.model.yrFinal   >= 2030 && xLim(2) < 2030, xLim(2) = 2050; end

      tscale = TspanScale(xLim,xLim);

      for i = 1:length(ax)
        fig.CurrentAxes = ax(i); pause(0.1)
        if numCol >= 5, ax(i).XTickLabelRotation = 45; end
        if strcmp(tscale.flag,'on') && xLim(1) <  tscale.Lim(2)
          tscale.SetTspanAxisScaled(xLim);
          if strcmp(obj.isUseSameAxis,"on") && i > 1
            ax(i).XTickLabel = ''; continue; 
          end
          xlim([tscale.tspan]);
        else
          xlim(xLim);
          xticks('auto');
        end
        if numCol >= 4 && xLim(1) <= 1750
          for j = 1:length(ax(i).XTickLabel)
            ax(i).XTickLabelRotation = 45;
            if length(ax(i).XTickLabel{j}) == 2,  ax(i).XTickLabel{j} = ''; end
          end
        end
      end
    end

  end %end methods
  
  methods(Static)
%% 
    function colorArray = color(plotName,iDir)
      switch plotName
        case {'_EDGARv5','EDGARv5'} 
          colorArray = "#005AFF";
        case {'_CEDS','CEDS'}
          colorArray = "#FF4B00";
        case {'_EDGARv6','EDGARv6'} 
          colorArray = "#03AF7A";
        case "Hist", colorArray = [0.6 0.6 0.6]; 
        case {"IMAGE-ssp119","SSP1-1.9"}, colorArray = [61 105 225]./255; %blue
        case {"IMAGE-ssp126","SSP1-2.6"}, colorArray = [34 139 34]./255; %green
        case {"MESSAGE-GLOBIOM-ssp245","SSP2-4.5"}, colorArray = [115,66,41]./255; %brown
        case {"AIM-ssp370","SSP3-7.0"}, colorArray = [189,183,107]./255; %gold
        case {"REMIND-MAGPIE-ssp534-over","SSP5-3.4os"}, colorArray = [186,85,211]./255; %pink
        case {"REMIND-MAGPIE-ssp585","SSP5-8.5"}, colorArray = [255,127,80]./255; %orange

        otherwise
          switch iDir
          case 1, colorArray = [0    0.4470    0.7410];
          case 2, colorArray = [0.8500    0.3250    0.0980];
          case 3, colorArray = [0.9290    0.6940    0.1250];
          case 4, colorArray = [0.4940    0.1840    0.5560];
          case 5, colorArray = [0.4660    0.6740    0.1880];
          case 6, colorArray = [0.3010    0.7450    0.9330];
          case 7, colorArray = [0.6350    0.0780    0.1840];
          end
      end
    end

    %%
    function setYAxisLimitCategoryTime(plotType,category,xLim,model)
      ax = gca; 

      switch plotType
      case {'ECH4','ECH4life'}
        switch category
          case 'tot'
            ax.YLim  = [200 750];
            ax.YTick = 200:50:750;
          case 'bio'
            ax.YLim  = [180 550];
            ax.YTick = 200:50:550;
          case {'anth_bio','natr_bio'}
            ax.YLim  = [0 370];
            ax.YTick = 0:50:370;      
          case 'ff'
            ax.YLim  = [0 250];
            ax.YTick = 0:40:250;
          case 'anth_ff'
            ax.YLim  = [0 210];
            ax.YTick = 0:40:200;
          case 'geo'
            ax.YLim  = [0 60];
            ax.YTick = 0:10:60;
          case 'bb'
            ax.YLim = [0 120];
            ax.YTick = 0:20:120;
          case 'rice'
            ax.YLim = [0 80];
            ax.YTick = 0:20:80;
          case 'rumi'
            ax.YLim = [0 180];
            ax.YTick = 0:20:180;
          case 'wast'
            ax.YLim = [0 120];
            ax.YTick = 0:20:120;
          case 'coal'
            ax.YLim = [0 120];
            ax.YTick = 0:20:120;
          case 'gas'
            ax.YLim = [0 120];
            ax.YTick = 0:20:120;
          case 'life'
            ax.YLim = [7.4 10.1];
            ax.YTick = 7.5:0.5:10;
        end
      case 'SourceFraction'
        switch category
          case 'bio'
            ax.YLim  = [60 100];
            ax.YTick = 60:5:100;
          case 'ff'
            ax.YLim  = [0 40];
            ax.YTick = 0:5:40;
          case 'anth_ff'
            ax.YLim  = [0 40];
            ax.YTick = 0:5:40;
          case 'geo'
            ax.YLim  = [0 20];
            ax.YTick = 0:5:20;
          case 'bb'
            ax.YLim  = [0 20];
            ax.YTick = 0:5:20;
        end
      case 'Parameter'
        [ymin,ymax]  = MyFunc.SetParaMinMax(category,model);
        ax.YLim      = [ymin ymax];
      case 'ParaTimeVar' 
          ax.YLim    = [0 model.paraIAVpercent.def];
      case 'Atmosphere'
        [ymin,ymax]  = MyFunc.SetAtmosMinMax(category,xLim(1),xLim(2)); 
        ax.YLim      = [ymin ymax];  
        if strcmp(category,'CH4') && xLim(1) < 1750
          ax.YTick = 800:400:2400; 
        elseif strcmp(category,'CH4')
          ax.YTick = 1500:100:1900;
          ax.YLim(2) = 1900;
        end 
      otherwise
        switch category
          case 'd13Ctot'
            ax.YLim = [-58 -51];
            ax.YTick = -58:1:-51;        
          case 'dDtot'
            ax.YLim = [-330 -250];
            ax.YTick = -330:10:-250;  
          case 'KIEC' 
            ax.YLim = [1.005 1.008];
            ax.YTick = 1.005:0.0005:1.008;   
            ax.YTickLabel = {'1.005','','1.006','','1.007','','1.008'};
          case 'KIED'
            ax.YLim = [1.25 1.30];
            ax.YTick = 1.25:0.01:1.30;         
            ax.YTickLabel = {'1.25','1.26','1.27','1.28','1.29','1.30'};  
        end
      end
    end

  end
end
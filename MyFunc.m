% A class that aggregates the functions for figure and text processing.
% Ryo Fujita 2024

classdef MyFunc 
  methods(Static)
    function [data,header] = readText(filename,colEnd,logicHeader)
      fid = fopen(filename, 'r');
      
      if fid == -1
        error('The file ' + filename + ' cannot be opened!');
      end
      
      delimiter = '-';
      line = fgetl(fid); 
      while ischar(line)
        if ~isempty(line) && strcmp(line(1),delimiter)
          break; 
        end
        line = fgetl(fid);
      end
      
      data = [];

      if logicHeader
        header = string(textscan(fid, repmat('%s', 1, colEnd),1));
      end

      if ischar(line) 
        data = textscan(fid, repmat('%f', 1, colEnd));
      end
      
      fclose(fid);
      
      if ~isempty(data)
        data = cell2mat(data);
      else
        error('No data is read!');
      end
    end

%% Opens datafile and reads time series data for free format and free data row 
%  Modified by Ryo Fujita in Dec 2018 from original code created by Heather Graven
    function varargout = readTextFreefmt(filename, format, delimiter)
      fid1=fopen(filename,'r');
      
      % Read strings delimited by a carriage return
      Header=textscan(fid1,'%s',1,'delimiter','\n'); 
      h=1;
      while h==1
          % Read strings delimited by a carriage return until '-' is encountered
          Header=textscan(fid1,'%s',1,'delimiter','\n'); 
          HeaderString=char(Header{1,1}(1));
          if numel(HeaderString)>0
              if strcmp(HeaderString(1),'-')==1 
              h=2;
              end
          end
      end
      InputText=textscan(fid1, format,'delimiter', delimiter); % Read data
      fclose(fid1);    
      
      if size(strsplit(format,' '),2) == size(InputText,2)     %strsplit
          for i = 1:1:size(InputText,2)
              varargout{i}=InputText{1,i};
          end
      else
          "Incorrect format setting!!"
          return
      end
    end

%% 
    function h = plotArea(tspan,y) 
      if isempty(y) , h = {}; return; end 

      if isa(y,'struct')
        categoryList = string(fieldnames(y)).';
      else
        categoryList = "None";
      end

      for catName = categoryList
        yPlot = getData(y,catName);
        if size(yPlot,2) == 1, h = plot(tspan,yPlot); return; end 
        [bottom,top] = bounds(yPlot,2);
        if length(tspan) == 1, h = []; return; end 
        h = area(tspan,[bottom top-bottom]); hold on
      set(h,'ShowBaseLine','off');
      set(h(1),'FaceColor','None','LineStyle','None');
      set(h(2),'EdgeAlpha',0.3,'LineStyle','-','LineWidth',1.3,'FaceAlpha',0.45);
      end

      function y = getData(y,catName)
        switch catName
        case 'None'
          ;
        otherwise
          y = y.(catName);
        end
      end
    end

%% 
    function subplotFront(obj,iCategory)
      if length(obj.categoryList) == 1, return; end
      subplot(obj.subPlotLayout(1),obj.subPlotLayout(2),iCategory);
    end    

%%
    function judge  = GetJudgeVaragin(obj,nameID)
      numRun  = size(obj,1);
      numNode = size(obj,2);
      judge(numRun,numNode).isMatch.all = 1;
      for i = 1:numRun
        for j = 1:numNode
          numCase                 = size(obj(i,j).(nameID),2);
          judge(i,j).isMatch.all  = true(1,numCase);
          judge(i,j).numMatch.all = sum(judge(i,j).isMatch.all);
        end
      end
    end

%% 
    function y_store = setPlotAtmosTarget(categoryList,subPlotLayout,model,xlim) 
      fig = gcf;
      y_store = []; %create array for store 
      for i = 1:length(categoryList)
        if ~isempty(subPlotLayout) && length(categoryList) > 2
          if contains(class(fig.Children),"TiledChartLayout")
            if i > 1, nexttile; end
          else
            subplot(subPlotLayout(1),subPlotLayout(2),i); 
          end
        else
          subplot(subPlotLayout(1),subPlotLayout(2),i)
        end

        category = categoryList(i);
        h = MyFunc.plotAtmosTarget(model,category,xlim); 
      end    

      axs = findobj(fig,'type','ax'); fig.CurrentAxes = axs(end);
    end

%% 
    function  h = plotAtmosTarget(model,gasID,xLim)
      for j = 1:2
        yrTarget    = zeros(1,model.numNode-2);
        targetMean  = zeros(1,model.numNode-2);
        targetRange = zeros(1,model.numNode-2);
        for iTerm = 1:model.numNode-2
          if model.atmosTargetList(iTerm).(gasID).flag == 0 ...
            || (model.atmosTargetList(iTerm).(gasID).min == -1000 ...
            && model.atmosTargetList(iTerm).(gasID).max == 1000)
            yrTarget(iTerm)    = NaN;
            targetMean(iTerm)  = NaN;
            targetRange(iTerm) = NaN;
          elseif (j == 1 && model.atmosTargetList(iTerm).tspan(2) < 1970) || ...
            (j == 2 && model.atmosTargetList(iTerm).tspan(2) >= 1970)
            yrTarget(iTerm)    = model.atmosTargetList(iTerm).tspan(2);
            atmosTarget        = model.atmosTargetList(iTerm);
            targetMean(iTerm)  = mean([atmosTarget.(gasID).min,atmosTarget.(gasID).max]);
            targetRange(iTerm) = abs(atmosTarget.(gasID).min-atmosTarget.(gasID).max)/2;
          end
        end
        
        if j == 1
          h = errorbar(TspanScale(yrTarget,xLim).tspan, targetMean, targetRange, ...
            'k.','Linewidth',0.7,'CapSize',4,'LineStyle','none','Marker','none'); 
          hold on
        else
          plot(TspanScale(yrTarget,xLim).tspan, targetMean+targetRange, ...
            'k_','Linewidth',0.7,'MarkerSize',3);
          plot(TspanScale(yrTarget,xLim).tspan, targetMean-targetRange, ...
            'k_','Linewidth',0.7,'MarkerSize',3);
        end
      end
    end

%%
    function setTextCategory(category)
      ax = gca;
      switch category
        case 'tot'
          txt = 'Total';
        case 'bio'
          txt = 'Biogenic';
        case 'anth_bio'
          txt = 'Anthropogenic Biogenic';
        case 'natr_bio'
          txt = 'Natural Biogenic';
        case 'ff'
          txt = 'Fossil Fuel';
        case 'anth_ff'
          txt = 'Anthropogenic Fossil Fuel';
        case 'geo'
          txt = 'Geologic';
        case 'bb'
          txt = 'Biomass Burning';
        case 'rice'
          txt = 'Rice';
        case 'rumi'
          txt = 'Ruminant';
        case 'wast'
          txt = 'Waste';
        case 'coal'
          txt = 'Coal';
        case 'gas'
          txt = 'Gas&Oil';
        case 'life'
          txt = 'Lifetime';
        case 'npr'
          txt = 'Nuclear ^{14}CH_{4}';
        otherwise 
         return   
      end 

      h_text = text(ax,0.03,0.95,txt,'Units','normalized','VerticalAlignment','top'); hold on %2021.6.3 add ax 2021.7.13 hold on
    end

%%
    function SetLabelTimeseries(namePlotList,propID)
      numName = size(namePlotList,2);
      if size(namePlotList,2) > 1, numSub = MyFunc.gridXYSubplot(namePlotList); end 

      for iName = 1:numName
        nameID = namePlotList(iName);
        if size(namePlotList,2) > 1, subplot(numSub(1),numSub(2),iName); end
        box on, hold on
        xlabel('Year')
        ylabel(MyFunc.nameLabel(nameID,propID))
      end
    end

%% 
    function SetLabelHistogram(namePlotList,propName,varargin)
      numName = size(namePlotList,2);
      numSub = MyFunc.gridXYSubplot(namePlotList); 

      if ~isempty(varargin) && strcmp(varargin{1},"tile")
        t = tiledlayout(numSub(1),numSub(2));
      end

      for iName = 1:numName
        nameID = namePlotList(iName);

        if ~isempty(varargin) && strcmp(varargin{1},"tile")
          nexttile
        else
          subplot(numSub(1),numSub(2),iName); 
        end

        box on, hold on
        title(MyFunc.nameLabel(nameID,propName))
      end
    end

%% 
    function y = nameLabel(nameIDList,propID,varargin)
      if ~isempty(varargin), optName = varargin{1}; else, optName = []; end

      y = cell(1,length(nameIDList));
      i = 0;
      for nameID = nameIDList 
        i = i + 1;
        switch nameID
        case 'CH4'
          switch propID
          case {'atmos','Atmosphere'}
            y(i) = {'CH_4 (ppb)'};
          case {'ems','Emission','ECH4'} 
            y(i) = {'E_{CH_4} (Tg/yr)'};
          case {'judge'} 
            y(i) = {strcat('N_{CH_4}')};
          otherwise
            y(i) = {'CH_4'};
          end  
        case 'd13C'
          switch propID
          case {'atmos','Atmosphere','ems','Emission','ECH4'} 
            y(i) = {strcat('\delta^{13}C-CH_4 (',char(8240),')')};
          case {'judge'}
             y(i) = {strcat('N_{\delta^{13}C-CH_4}')};
          otherwise
            y(i) = {'\delta^{13}C-CH_4'};
          end    
        case 'dD'
          switch propID
          case {'atmos','Atmosphere','ems','Emission','ECH4'} 
             y(i) = {strcat('\deltaD-CH_4 (',char(8240),')')};
          case {'judge'}
             y(i) = {strcat('N_{\deltaD-CH_4}')};
          otherwise
            y(i) = {'\deltaD-CH_4'};
          end  
        case 'D14C'
          switch propID
          case {'atmos','Atmosphere','ems','Emission','ECH4'} 
             y(i) = {strcat('\Delta^{14}C-CH_4 (',char(8240),')')};
          case {'judge'}
             y(i) = {strcat('N_{\Delta^{14}C-CH_4}')};
          otherwise
            y(i) = {'\Delta^{14}C-CH_4'};
          end
        case 'dD14Cdt' 
          y(i) = {strcat('\Delta^{14}C-CH_4 growth (',char(8240),' yr^{-1})')};
        case 'all'
          switch propID
          case {'judge'}
             y(i) = {strcat('N_{All}')};
          end
        case 'C13'
          switch propID
          case {'atmos','Atmosphere'}
            y(i) = {'^{13}CH_4 (ppb)'};
          case {'ems','Emission','ECH4'} 
            y(i) = {'E_{C13} (Tg/yr)'};
          end
        case 'CHD'
          switch propID
          case {'atmos','Atmosphere'}
            y(i) = {'CH_3D (ppb)'};
          case {'ems','Emission','ECH4'} 
            y(i) = {'E_{CHD} (Tg/yr)'};
          end
        case 'C14'
          switch propID
          case {'atmos','Atmosphere'}
            y(i) = {'^{14}CH_4 (ppb)'};
          case {'ems','Emission','ECH4'}
            y(i) = {'E_{C14} (Tg/yr)'};
          end
        case 'fbb'
            y(i) = {'f_{bb}'};
        case 'fbio'    
            y(i) = {'f_{bio}'};         
        case 'fff'    
            y(i) = {'f_{ff}'};
        case 'Egeo'
            y(i) = {'E_{geo} (TgCH_{4}yr^{-1})'};
        case 'fanth_ff'
            y(i) = {'f_{anth\_ff}'};
        case 'fanth_bio'
            y(i) = {'f_{anth\_bio}'};
        case 'fnatr_bio'
            y(i) = {'f_{natr\_bio}'};
        case 'd13Ctot'
            y(i) = {strcat('\delta^{13}C_{tot} (',char(8240),')')};
        case 'd13Cbb'
            y(i) = {strcat('\delta^{13}C_{bb} (',char(8240),')')};
        case 'd13Cbio'
            y(i) = {strcat('\delta^{13}C_{bio} (',char(8240),')')};
        case 'd13Canth_bio'
            y(i) = {strcat('\delta^{13}C_{anth\_bio} (',char(8240),')')}; 
        case 'd13Cnatr_bio'
            y(i) = {strcat('\delta^{13}C_{natr\_bio} (',char(8240),')')};             
        case 'd13Canth_ff'
            y(i) = {strcat('\delta^{13}C_{anth\_ff} (',char(8240),')')};
        case 'd13Cgeo'
            y(i) = {strcat('\delta^{13}C_{geo} (',char(8240),')')};
        case 'd13Crumi'
            y(i) = {strcat('\delta^{13}C_{rumi} (',char(8240),')')};
        case 'd13Cff'
            y(i) = {strcat('\delta^{13}C_{ff} (',char(8240),')')};
        case 'd13CO2'
            y(i) = {strcat('Atmospheric \delta^{13}CO_{2} (',char(8240),')')};
        case 'dDtot'
            y(i) = {strcat('\deltaD_{tot} (',char(8240),')')};
        case 'dDbb'
            y(i) = {strcat('\deltaD_{bb} (',char(8240),')')};
        case 'dDbio'
            y(i) = {strcat('\deltaD_{bio} (',char(8240),')')};
        case 'dDanth_bio'
            y(i) = {strcat('\deltaD_{anth\_bio} (',char(8240),')')}; 
        case 'dDnatr_bio'
            y(i) = {strcat('\deltaD_{natr\_bio} (',char(8240),')')}; 
          case 'dDanth_ff'
            y(i) = {strcat('\deltaD_{anth\_ff} (',char(8240),')')};
        case 'dDgeo'
            y(i) = {strcat('\deltaD_{geo} (',char(8240),')')};
        case 'dDff'
            y(i) = {strcat('\deltaD_{ff} (',char(8240),')')};
        case 'tau'
            y(i) = {'\tau_{bios} (yr)'};
        case 'phi'
            y(i) = {'\Phi (GBq/GWa)'};
        case 'life'
            y(i) = {'\lambda_{tot} (yr)'};
        case 'KIEC'
            y(i) = {'KIE^{C}'};
        case 'KIED'
            y(i) = {'KIE^{D}'};
        case {'tot','bio','ff','geo','bb'}
          switch propID
          case {'SourceWorking','SourceSecondary','SourcePrimary','ems',...
            'Emission','ECH4','ECH4Life','ECH4NPP'} 
            y(i) = {strcat('E_{',char(nameID),'} (Tg/yr)')};
          case 'SourceFraction'
            y(i) = {strcat('Frac_{',char(nameID),'} (%)')};
          otherwise
            y(i) = {char(nameID)};
          end
        case 'npr'
          switch propID
            case 'ECH4NPP'
              y(i) = {'Nuclear ^{14}CH_{4} (TBq/yr)'};
          end
        case 'anth_ff'
          switch propID
          case {'SourceWorking','SourceSecondary','SourcePrimary','ems','Emission','ECH4','ECH4Life'} 
            y(i) = {strcat('E_{anth\_ff} (Tg/yr)')};
          case 'SourceFraction'
            y(i) = {'Frac_{anth\_ff} (%)'};
          otherwise
            y(i) = {'anth\_ff'};
          end
        case 'anth_bio'
          switch propID
          case {'SourceWorking','SourceSecondary','SourcePrimary','ems','Emission','ECH4','ECH4Life'} 
            y(i) = {strcat('E_{anth\_bio} (Tg/yr)')};
          case 'SourceFraction'
            y(i) = {'Frac_{anth\_bio} (%)'};
          otherwise
            y(i) = {'anth\_bio'};
          end
        case 'natr_bio'
          switch propID
          case {'SourceWorking','SourceSecondary','SourcePrimary','ems','Emission','ECH4','ECH4Life'} 
            y(i) = {strcat('E_{natr\_bio} (Tg/yr)')};
          case 'SourceFraction'
            y(i) = {'Frac_{natr\_bio} (%)'};
          otherwise
            y(i) = {'natr\_bio'};
          end       
        case 'tspan'
          y(i) = {'Year'};
        otherwise
            y(i) = {char(nameID)};
        end
      end

      if ~isempty(optName) 
        switch optName
        case 'cutUnit'
          y = MyFunc.cutUnit(y);
        otherwise
          error('Proper option name for nameLabel is not set')
        end
      end
    end

%% 
    function y = cutUnit(labelList)
      y = labelList;
      i = 0;
      for label = labelList
        i = i + 1;
        if contains(labelList(i)," ")
          y(i) = extractBefore(labelList(i)," ");
        end
      end
    end

%%
    function y = labelAnomaly(AnomSet,y_org)
      y = strcat('Anomaly of'," ",y_org);
      
      if strcmp(AnomSet.type,'percent')
        y = extractBefore(y,'(');
        y = cellstr(strcat(y,' (%)'));
      end
    end

%%
    function SetAxisLimitTimeseries(namePlotList,propID,model)
      numName  = size(namePlotList,2);
      if size(namePlotList,2) > 1, numSub = MyFunc.gridXYSubplot(namePlotList); end 
      for iName = 1:numName
        nameID = namePlotList(iName);
        if size(namePlotList,2) > 1, subplot(numSub(1),numSub(2),iName); end
        ax = gca; 
        
        ax = MyFunc.SetAxisLimit(ax,"YAxis",nameID,propID,model);
        ax = MyFunc.SetXAxisTime(model,ax); 
        ax.XMinorTick = 'on';ax.YMinorTick = 'on';
      end
    end

%%
    function ax = SetXAxisTime(model,ax)
      if model.yrInitial == 1750
        xmin = 1680;
      else
        xmin = floor((model.yrInitial - model.timeSpinUpYr)/10 ) * 10; 
      end

      xmax = ceil(model.yrFinal/10 ) * 10; 

      ax.XAxis.Limits = [xmin xmax]; 
    end

%% 
    function AX = SetAxisLimit(AX,nameAxis,nameID,propID,model)
      switch propID
      case {'atmos','Atmosphere'}
        [ymin, ymax] = MyFunc.SetAtmosMinMax(nameID,model.yrInitial,model.yrFinal); 
      case {'ems','Emission'}
        [ymin, ymax] = MyFunc.SetEmsMinMax(nameID); 
      case {'para','ParameterSourceSink'}  
        [ymin, ymax] = MyFunc.SetParaMinMax(nameID,model);
      case {'judge','JudgeAtmosToAtmosTarget'}  
        ymin     = 0;
        ymax     = model.numCase_org;
      case {'SourceWorking','SourceSecondary','SourcePrimary'} 
        [ymin, ymax] = MyFunc.SetSourceMinMax(nameID);
      case {'SourceFraction'}
        [ymin, ymax] = MyFunc.SetSourceFractionMinMax(nameID); 
      otherwise
        return
      end

      AX.(nameAxis).Limits = [ymin ymax]; 
    end

%% 
    function y = EdgePara(plotClass,model,nameID)
      switch plotClass
      case {'atmos','Atmosphere'}
        [ymin, ymax] = MyFunc.SetAtmosMinMax(nameID,model.yrInitial,model.yrFinal);
      case {'ems','Emission'}
        [ymin, ymax] = MyFunc.SetEmsMinMax(nameID);          
      case {'para','ParameterSourceSink'}  
        [ymin, ymax] = MyFunc.SetParaMinMax(nameID,model); 
      case {'SourceWorking','SourceSecondary','SourcePrimary'}
        [ymin, ymax] = MyFunc.SetSourceMinMax(nameID); 
      case {'SourceFraction'}
        [ymin, ymax] = MyFunc.SetSourceFractionMinMax(nameID); 
      case 'ParaTimeVar'
        ymin = 0;
        ymax = model.paraIAVpercent.def;
      otherwise 
        ax = gca;
        ymin = ax.XLim(1);
        ymax = ax.XLim(2);
      end

      di = (ymax - ymin)/20;
      y  =  ymin  : di : ymax;
    end

%%
    function [ymin, ymax] = SetAtmosMinMax(nameID,yrInitial,yrFinal)
        switch nameID
        case 'CH4' 
          if yrFinal <= 2025
            ymax = 2400; 
          else
            ymax = 3500;
          end  

          if yrInitial <= 1750 
            ymin = 600; 
          elseif yrInitial <= 1850
            ymin = 700; 
          elseif yrInitial <= 1900
            ymin = 800; 
          elseif yrInitial <= 1950
            ymin = 1000; 
          elseif yrInitial <= 1960
            ymin = 1200; 
          elseif yrInitial <= 1970
            ymin = 1300; 
          elseif yrInitial <= 1980
            ymin = 1500;
            ymax = 1900; 
          else
            ymin = 1600;
          end

        case 'd13C'
          if yrFinal <= 2025
            ymax = -45.1;
          else
            ymax = -45;
          end

          if yrInitial <= 1950 
            ymin = -52.9;
          elseif yrInitial <= 1970
            ymin = -50; 
          else
            ymin = -48.5;
            ymax = -46;
          end

        case 'dD'
          if yrInitial <= 1750 
            ymin = -130; 
          else
            ymin = -120;
          end

          if yrFinal <= 2025
            ymax = -60;
          else
            ymax = -20; 
          end
        case 'D14C'
          if yrInitial <= 1950 
            ymin = -300; 
          elseif yrInitial <= 1970 
            ymin = -200; 
          elseif yrInitial <= 1980 
            ymin = 0; 
          else
            ymin = 100; 
          end

          if yrFinal <= 2025
            ymax = 500;
          else
            ymax = 1000; 
          end
        end           
    end

%%
    function [ymin, ymax] = SetEmsMinMax(nameID)
        switch nameID
        case 'CH4'
          ymin = 200; ymax = 800;
        case {'d13C','d13Ctot'}
          ymin = -65; ymax = -50;
        case {'dD','dDtot'}
          ymin = -320; ymax = -250;
        case 'D14C'
          ymin = -300; ymax = 800;
        end 
    end

%% 
    function [ymin, ymax] = SetSourceMinMax(nameID)
      switch nameID
        case 'tot'
          ymin = 180; ymax = 700;
        case 'bio'
          ymin = 180; ymax = 520;
        case 'anth_bio'
          ymin = 0; ymax = 310;           
        case 'natr_bio'
          ymin = 90; ymax = 310;           
        case 'ff'
          ymin = 0; ymax = 250;
        case 'anth_ff'
          ymin = 0; ymax = 210;
        case 'geo'
          ymin = 0; ymax = 60;
        case 'bb'
          ymin = 0; ymax = 60;
      end
    end

%% 
    function [ymin, ymax] = SetSourceFractionMinMax(nameID)
      switch nameID
        case 'bio'
          ymin = 50; ymax = 100;
        case 'anth_bio'
          ymin = 0; ymax = 60;
        case 'natr_bio'
          ymin = 30; ymax = 90;          
        case 'ff'
          ymin = 0; ymax = 40;
        case 'anth_ff'
          ymin = 0; ymax = 35;
        case 'geo'
          ymin = 0; ymax = 35;
        case 'bb'
          ymin = 0; ymax = 25;
      end
    end

%% 
    function [ymin, ymax] = SetParaMinMax(nameID,model)
      switch nameID
        case 'fbb'
          ymin = 0.5; ymax = 3.5;
        case 'fbio'
          ymin = 0.5; ymax = 1.5;
        case 'fanth_bio'
          ymin = 0.5; ymax = 1.5;
        case 'fnatr_bio'
          ymin = 0.5; ymax = 1.5;
        case 'fanth_ff'
          ymin = 0.5; ymax = 1.5;
        case 'Egeo'
          ymin = 0; ymax = 80;
        case 'tau'
          ymin = 1; ymax = 12;
        case 'phi'
          ymin = 80; ymax = 380;
        case 'life'
          ymin = 7.5; ymax = 10;
        case 'floss'
          ymin = 0.9; ymax = 1.1;
        case 'KIEC'
          ymin = 1.005; ymax = 1.008;
        case 'KIED'
          ymin = 1.25; ymax = 1.30; 
        case 'd13Cbio'
          ymin = -65.4; ymax = -59.0;
        case 'd13Canth_bio'
          ymin = -65.4; ymax = -59.0;
        case 'd13Cnatr_bio'
          ymin = -65.4; ymax = -59.0;
        case 'd13Cbb'
          ymin = -29.8; ymax = -14.6;
        case 'd13Cff'
          ymin = -46.8; ymax = -41.2;
        case 'd13Canth_ff'
          ymin = -46.8; ymax = -41.2;
        case 'd13Cgeo'
          ymin = -52.0; ymax = -46.0;
        case 'dDbio'
          ymin = -332; ymax = -302;
        case 'dDanth_bio'
          ymin = -332; ymax = -302;
        case 'dDnatr_bio'
          ymin = -332; ymax = -302;
        case 'dDbb'
          ymin = -226; ymax = -196;
        case 'dDff'
          ymin = -212; ymax = -182;
        case 'dDanth_ff'
          ymin = -212; ymax = -182;
        case 'dDgeo'
          ymin = -212; ymax = -182;
      end
    end

%% 
    function y = ColorCategory(cname)
      switch cname
        case 'bio'
          y = [0.2    0.4    1];
        case 'ff'
          y = [0.800    0.1250    0.2980];
        case 'anth_ff'
          y = [0.800    0.1250    0.1080];
        case 'geo'
          y = [1    0.8    0];
        case 'bb'
          y = [0    0.45    0];
        case 'npr'
          y = [0.5   0    0.5];
        case 'tot'
          y = [0   0    0];
      end  
    end

%% 
    function y = gridXYSubplot(namePlotList)  
      numName = length(namePlotList); 
      if numName == 0
        disp('No variable parameter')
        numRow = 0; numCol = 0;
      elseif numName == 1 
        numRow = 1; numCol = 1;
      elseif numName == 2 
        numRow = 2; numCol = 1;
      elseif 2 < numName && numName <= 4
        numRow = 2; numCol = 2;
      elseif 4 < numName && numName <= 6 
        numRow = 2; numCol = 3;
      elseif 6 < numName && numName <= 9 
        if numName == 8 && sum(ismember(namePlotList,'d13Cbio'),2) == 1
          numRow = 2; numCol = 4;
        else
          numRow = 3; numCol = 3;
        end
      elseif 9 < numName && numName  <= 12 
        if numName == 10  && sum(ismember(namePlotList,'d13Canth_bio'),2) == 1
          numRow = 2; numCol = 5;
        else
          numRow = 3; numCol = 4;
        end        
      elseif 12 < numName && numName  <= 16 
        numRow = 4;
        numCol = 4;
      elseif 16 < numName && numName  <= 20 
        numRow = 4;
        numCol = 5;
      elseif 20 < numName && numName  <= 25 
        numRow = 5;
        numCol = 5;
      end
      y = [numRow numCol];
    end

%%
    function y = ColorNode(iNode)
      if iNode == 1; y = [0.5 0.5 0.5]; end                 %dark gray
      if iNode == 2; y = [0    0.4470    0.7410]; end       %blue
      if iNode == 3; y = [0.8500    0.3250    0.0980]; end  %yeliow
      if iNode == 4; y = [0.9290    0.6940    0.1250]; end  %orange
      if iNode == 5; y = [0.4940    0.1840    0.5560]; end
    end

%%
    function c = hsv_v2(numNode)
      c = hsv(numNode);
      if numNode == 11
        c(3,:)  = [0.82,0.82,0.08];
        c(5,:)  = [0,0.5,0];
        c(6,:)  = [0.125,0.695,0.66];
        c(11,:) = [0.777,0.08,0.519];
      elseif numNode == 12 %2020.12.3 
        c(3,:)  = [0.82,0.82,0.08];
        c(5,:)  = [0,0.5,0];
        c(6,:)  = [0.125,0.695,0.66];
        c(11,:) = [0.9098 	0.3216 	0.5961];
        c(12,:) = [0.777,0.08,0.519 ]; 
      end  
    end
  end

end



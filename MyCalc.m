% A class that aggregates the functions for data processing.
% Ryo Fujita 2024

classdef MyCalc 
  methods(Static)
    function plotObj = SetPlotObj(plotType,model,varargin)
      if ~isempty(varargin), nameList = varargin{1}; else, nameList = []; end
      switch plotType 
        case 'Atmosphere'
          if isempty(model.atmos)
            plotObj = Atmosphere(model); 
          else
            plotObj = model.atmos; 
          end
          
          if model.status ~= 10 
            plotObj.CH4 = plotObj.CH4 ./ MassConstants.factor; 
          end
        case 'Parameter'
          if isempty(model.para)
            plotObj = ParameterSourceSink(model);
          else
            plotObj = model.para; 
          end
          if ~isempty(nameList), plotObj = MyCalc.AddCalc(plotObj,nameList,model); end 
        otherwise
          if isempty(model.para), model.para = ParameterSourceSink(model); end 
          model.numCase = size(model.para.fbb,2); 
          model.numData = length(model.para.tspan);
        switch plotType  
          case "ECH4"
            plotObj = Emission().getCH4only(model,nameList); 
          case "ECH4NPP" 
            M       = MassConstants;
            plotObj = Emission().getCH4only(model,nameList); 
            plotObj.npr = Emission().SetEmsINodeRadiocarbon(model,model.para,"npr").C14.npr; 
            plotObj.npr = plotObj.npr./M.W_CH4./M.Bq;%Tg(CH4)yr^-1 -> Tmole(14C)yr^-1 -> TBq(14CH4)yr^-1
          case "SourceFraction"
            plotObj = SourceFraction(model,nameList); 
          otherwise
            if isempty(nameList), error('Need to set nameList to use Struct class!'); end
            plotObj.tspan = model.para.tspan; 
            for name = nameList
              switch name
              case {'d13Ctot','dDtot'} 
                E_CH4 = Emission().getCH4only(model,["anth_bio","natr_bio","anth_ff","geo","bb"]);
                plotObj.(name) = SourceWorking.(name)(model.para,E_CH4); 
              otherwise 
                if ismember(name,model.nameList.para)            
                  plotObj.(name) = model.para.(name);              
                elseif ismember(name,string(fieldnames(SourceWorking)).')
                  plotObj.(name) = Emission().getCH4only(model,name).(name);              
                elseif ismember(name,["d13Cbio","d13Cff","dDbio","dDff"])
                  plotObj.(name) = MyCalc.getDeltaBIOorFF(name,model); 
                elseif strcmp(name,"life") && ~isempty(model.para.floss)
                  plotObj.(name) = MyCalc.convertFLossToLife(model); 
                elseif strcmp(name,"ohAnomaly") && ~isempty(model.para.floss) 
                  plotObj.(name) = MyCalc.convertFLossToOHAnomaly(model); 
                elseif strcmp(name,"totLoss") && ~isempty(model.para.floss)
                  plotObj.(name) = MyCalc.convertFLossToLoss(model); 
                elseif strcmp(name,"oh") && ~isempty(model.para.floss)
                  plotObj.(name) = MyCalc.convertFLossToOH(model); 
                elseif strcmp(name,"ohLife") && ~isempty(model.para.floss)
                  plotObj.(name) = 1./MyCalc.convertFLossToOH(model); 
                else
                  error('No such a property name both in SourceWorking or ParameterSourceSink!');
                end
              end
            end
        end
      end
    end

%% 
    function y = getBinLimitsList(plotType,calcType,nameList,model,obj2D)
      for nameID = nameList
        switch plotType       
        case {'Parameter'}  
          if ismember(nameID,model.nameList.para) && ~strcmp(calcType,"PeriodDiff")
            minLim = model.paraTargetList(end).(nameID).min;
            maxLim = model.paraTargetList(end).(nameID).max;
          else
            [minLim, maxLim] = bounds(obj2D,'all');
          end
        otherwise 
          [minLim, maxLim] = bounds(obj2D,'all');
        end
        y.(nameID) = [minLim, maxLim];
      end
    end

%%
    function y = binLimits(obj2D)
      if size(obj2D,1) == 1
        y = obj2D; return
      end

      [minLim, maxLim] = bounds(obj2D,'all');
      y = [minLim, maxLim];
    end

%% 
    function y = values(h,h_new)
      if length(h) == 1
        y = h.Values;
        return
      end

      switch class(h)
        case 'matlab.graphics.chart.primitive.Histogram'
          set(h,'BinEdges',h_new.BinEdges);
          y = mean(vertcat(h.Values),1);
        case 'matlab.graphics.chart.primitive.Histogram2'
          set(h,'XBinEdges',h_new.XBinEdges);
          set(h,'YBinEdges',h_new.YBinEdges);
          y_tmp = reshape(horzcat(h.Values),h_new.NumBins(1),h_new.NumBins(2),length(h)); 
          y     = mean(y_tmp,3);
      end
    end

%% 
    function y = binEdges(lim,numBins)
      if lim(2) == lim(1)
        y    = zeros(1,numBins+1);
        y(:) = lim(1);
      else
        res = (lim(2) - lim(1))/numBins;
        y   = lim(1):res:lim(2);
      end
    end

%%
    function y = meanPropability(obj2D,binEdges)
      prop = zeros(length(binEdges)-1,size(obj2D,2));

      for i = 1:size(prop,2)
        prop(:,i) = histcounts(obj2D(:,i),binEdges, ...
        'Normalization', 'probability');
      end
      y = mean(prop,2);
    end

%% 
    function  y = getStatsHistogram(y,varargin)
      if isa(y,'double') && size(y,2) == 2
          data = y;
          y = varargin{1};
          y = getStatsData2(data,y);
      elseif y.plotStyleType == "Histogram"
          y = getStatsHist(y);
      elseif y.plotStyleType == "Histogram2"
          y = getStatsHist2(y);
      else
        error("Not yet developed!")
      end

      function y = getStatsData2(data,y)
        y.Xave = mean(data(:,1));  
        y.Xstdev = std(data(:,1)); 
        y.Xcdf = [];  
        y.Xmed = median(data(:,1),1);
        y.Xcdf50 = 0.5; 
        y.Xprc68 = prctile(data(:,1),[16,84]);
        y.Xcdf68 = [0.16,0.84];
        y.Xprc95 = prctile(data(:,1),[2.5,97.5]);
        y.Xcdf95 = [0.025,0.975];

        y.Yave = mean(data(:,2));     
        y.Ystdev = std(data(:,2));  
        y.Ycdf = [];   
        y.Ymed = median(data(:,2),1);
        y.Ycdf50 = 0.5; 
        y.Yprc68 = prctile(data(:,2),[16,84]);
        y.Ycdf68 = [0.16,0.84];
        y.Yprc95 = prctile(data(:,2),[2.5,97.5]);
        y.Ycdf95 = [0.025,0.975];
 
        R = corrcoef(data(:,1),data(:,2));
        y.R = R(1,2);
      end

      function y = getStatsHist(y)
        assert(length(y) == 1, 'Need to input single object!');
        assert(~isempty(y.Values)); 
        y.ave   = MyCalc.average(y);
        y.stdev = MyCalc.stdev(y,y.ave); 
        y.cdf   = MyCalc.getCdf(y); %figure; histogram('BinEdges',h.BinEdges,'BinCounts',h.Values); figure; plot(h.BinEdges,y.cdf);
        [y.med,  y.cdf50] = MyCalc.med(y);
        [y.prc68,y.cdf68] = MyCalc.prc68(y);
        [y.prc95,y.cdf95] = MyCalc.prc95(y);
      end

      function y = getStatsHist2(y)
        assert(length(y) == 1, 'Need to input single object!');
        assert(~isempty(y.Values)); 
        y.Xave   = MyCalc.average(y,'X');
        y.Xstdev = MyCalc.stdev(y,y.Xave,'X'); 
        y.Xcdf   = MyCalc.getCdf(y,'X');
        [y.Xmed,  y.Xcdf50] = MyCalc.med(y,'X');
        [y.Xprc68,y.Xcdf68] = MyCalc.prc68(y,'X');
        [y.Xprc95,y.Xcdf95] = MyCalc.prc95(y,'X');

        y.Yave   = MyCalc.average(y,'Y');
        y.Ystdev = MyCalc.stdev(y,y.Yave,'Y');
        y.Ycdf   = MyCalc.getCdf(y,'Y');
        [y.Ymed,  y.Ycdf50] = MyCalc.med(y,'Y');
        [y.Yprc68,y.Ycdf68] = MyCalc.prc68(y,'Y');
        [y.Yprc95,y.Ycdf95] = MyCalc.prc95(y,'Y');
 
        y.R = MyCalc.RCorr(y); 
      end
    end

%%
    function y = getCdf(h,XorY)
      switch class(h)
        case 'matlab.graphics.chart.primitive.Histogram'
          y = getCdfHist(h);
        case 'matlab.graphics.chart.primitive.Histogram2'
          y = getCdfHist2(h,XorY);
        case 'struct' 
        switch h.plotStyleType
         case 'Histogram'
          y = getCdfHist(h);
         case 'Histogram2'
          y = getCdfHist2(h,XorY);
        end
      end
      
      function y = getCdfHist(h)
        y = [0 cumsum(h.Values)];
      end
      function y = getCdfHist2(h,XorY)
        switch XorY
         case 'X', y = [0 cumsum(sum(h.Values,2),1).'];
         case 'Y', y = [0 cumsum(sum(h.Values,1),2)];
        end
      end
    end

%%
    function y = average(h,XorY)
      switch class(h)
        case 'matlab.graphics.chart.primitive.Histogram'
          y = sum(movmean(h.BinEdges,2,'Endpoints','discard').*h.Values);
        case 'matlab.graphics.chart.primitive.Histogram2'
          switch XorY
            case 'X'
              y = sum(movmean(h.XBinEdges,2,'Endpoints','discard').*h.Values.','all');
            case 'Y'
              y = sum(movmean(h.YBinEdges,2,'Endpoints','discard').*h.Values,'all');
          end
        case 'struct'
          switch h.plotStyleType
            case 'Histogram'
              y = sum(movmean(h.BinEdges,2,'Endpoints','discard').*h.Values);
            case 'Histogram2'
            switch XorY
            case 'X'
              y = sum(movmean(h.XBinEdges,2,'Endpoints','discard').*h.Values.','all');
            case 'Y'
              y = sum(movmean(h.YBinEdges,2,'Endpoints','discard').*h.Values,'all');
            end
          end
      end
    end

%%
    function y = stdev(h,ave,XorY) 
      switch class(h)
        case 'matlab.graphics.chart.primitive.Histogram'
          y = stdevHist(h,ave);
        case 'matlab.graphics.chart.primitive.Histogram2'
          y = stdevHist2(h,ave,XorY);
        case 'struct'
         switch h.plotStyleType
         case 'Histogram'
          y = stdevHist(h,ave);
         case 'Histogram2'
          y = stdevHist2(h,ave,XorY);
        end
      end

      function y = stdevHist(h,ave)
        y = sum((movmean(h.BinEdges,2,'Endpoints','discard')-ave).^2.*h.Values);
      end

      function y = stdevHist2(h,ave,XorY)
       switch XorY
        case 'X', y = sqrt(sum((movmean(h.XBinEdges,2,'Endpoints','discard')-ave).^2.*h.Values.','all'));
        case 'Y', y = sqrt(sum((movmean(h.YBinEdges,2,'Endpoints','discard')-ave).^2.*h.Values,'all'));
        end
      end
    end

%%
    function [prc,frac] = prc50(h,varargin)
      if isempty(varargin)
        [prc,frac] = MyCalc.prcXX(h,0.50);
      else  
        [prc,frac] = MyCalc.prcXX(h,0.50,varargin{1});
      end
    end

%%
    function [prc,frac] = prc68(h,varargin)
      if isempty(varargin)
        [prc,frac] = MyCalc.prcXX(h,0.68);
      else  
        [prc,frac] = MyCalc.prcXX(h,0.68,varargin{1});
      end
    end

%%
    function [prc,frac] = prc90(h,varargin)
      if isempty(varargin)
        [prc,frac] = MyCalc.prcXX(h,0.90);
      else  
        [prc,frac] = MyCalc.prcXX(h,0.90,varargin{1});
      end
    end

%%
    function [prc,frac] = prc95(h,varargin)
      if isempty(varargin)
        [prc,frac] = MyCalc.prcXX(h,0.95);
      else  
        [prc,frac] = MyCalc.prcXX(h,0.95,varargin{1});
      end
    end

%% 
    function [prc,frac] = prcXX(h,prc,XorY)
      prc1 = 0.5 - prc/2;  
      prc2 = 0.5 + prc/2;
      switch class(h)
        case 'matlab.graphics.chart.primitive.Histogram'
          [prc,frac] = prcXXHist(h,prc1,prc2);
        case 'matlab.graphics.chart.primitive.Histogram2'
          [prc,frac] = prcXXHist2(h,prc1,prc2,XorY); 
        case 'struct' 
        switch h.plotStyleType
         case 'Histogram'
          [prc,frac] = prcXXHist(h,prc1,prc2);
         case 'Histogram2'
          [prc,frac] = prcXXHist2(h,prc1,prc2,XorY);
        end
        case 'double' 
          prc  = prctile(h,[prc1*100 prc2*100],2);
          frac = [];
      end
      
      function [prc,frac] = prcXXHist(h,prc1,prc2)
        if isa(h,'struct') && ~isempty(h.cdf), cdf = h.cdf; else, cdf = MyCalc.getCdf(h); end
        [~,idx1] = min(abs(cdf-prc1));
        [~,idx2] = min(abs(cdf-prc2));
        prc  = [h.BinEdges(idx1) h.BinEdges(idx2)];
        frac = [cdf(idx1) cdf(idx2)];
      end

      function [prc,frac] = prcXXHist2(h,prc1,prc2,XorY) 
        XYcdf = strcat(XorY,'cdf');
        if isa(h,'struct') && ~isempty(h.(XYcdf)), cdf = h.(XYcdf); else, cdf = MyCalc.getCdf(h,XorY); end

        [~,idx1] = min(abs(cdf-prc1));
        [~,idx2] = min(abs(cdf-prc2));
        frac = [cdf(idx1) cdf(idx2)];
        switch XorY
          case 'X'
            prc = [h.XBinEdges(idx1) h.XBinEdges(idx2)];
          case 'Y'
            prc = [h.YBinEdges(idx1) h.YBinEdges(idx2)];
        end
      end
    end

%%
    function [prc,frac] = med(h,XorY)
      switch class(h)
        case 'matlab.graphics.chart.primitive.Histogram'
          [prc,frac] = medHist(h);
        case 'matlab.graphics.chart.primitive.Histogram2'
          [prc,frac] = medHist2(h,XorY);
        case 'struct'
        switch h.plotStyleType
         case 'Histogram'
          [prc,frac] = medHist(h);
         case 'Histogram2'
          [prc,frac] = medHist2(h,XorY);
        end
      end
      
      function [prc, frac] = medHist(h)
        if isa(h,'struct') && ~isempty(h.cdf)
          cdf = h.cdf; 
        else
          cdf = MyCalc.getCdf(h); 
        end

        [~,idx1] = min(abs(cdf-0.5));
        prc  = h.BinEdges(idx1);
        frac = cdf(idx1);
      end

      function [prc, frac] = medHist2(h,XorY)
        XYcdf = strcat(XorY,'cdf');
        if isa(h,'struct') && ~isempty(h.(XYcdf))
          cdf = h.(XYcdf); 
        else
          cdf = MyCalc.getCdf(h,XorY); 
        end

        [~,idx1] = min(abs(cdf-0.5));
        switch XorY
         case 'X', prc = h.XBinEdges(idx1);
         case 'Y', prc = h.YBinEdges(idx1);
        end
        frac = cdf(idx1);
      end
    end

%% 
    function y = meanHist(h)
      y = mean(h.Data);
    end

%%
    function y = prctl95Hist(h)
      y(1) = prctile(h.Data,2.5);
      y(2) = prctile(h.Data,97.5);
    end

%%
    function y = prctl68Hist(h)
      y(1) = prctile(h.Data,16);
      y(2) = prctile(h.Data,84);
    end

%%
    function y = RCorr(h)
      switch class(h)
        case 'matlab.graphics.chart.primitive.Histogram2'
          [R,p] = corrcoef(h.Data(:,1),h.Data(:,2));
          y = R(1,2);  
        case 'struct'
          switch h.plotStyleType
          case 'Histogram2'
             cov_tmp = ((movmean(h.XBinEdges,2,'Endpoints','discard')-h.Xave)./h.Xstdev).' * ...
                       ((movmean(h.YBinEdges,2,'Endpoints','discard')-h.Yave)./h.Ystdev);
             y = sum(cov_tmp.*h.Values,'all'); 
          end
      end
    end

%% 
    function y = wtmeanTimeseries(obj1,obj2,varargin)
      if ~isequal(obj1.tspan,obj2.tspan), error('Need to set same tspan object'); end
      
      if ~isempty(varargin) 
        nameList = varargin{1};
      else
        nameList  = string(fieldnames(obj1)).';
        nameIndex = string(fieldnames(obj1)).' ~= "tspan";        
        nameList  = nameList(nameIndex);
      end

      y.tspan = obj1.tspan;
      numData = length(y.tspan);

      for name = nameList
        if isempty(obj1.(name)), continue; end
        y.(name).ave   = zeros(numData,1);
        y.(name).sd    = zeros(numData,1);
        y.(name).prc68 = zeros(numData,2);
        for i = 1:numData
          [y.(name).ave(i),y.(name).sd(i)] = MyCalc.wtMean_SD(...
            mean(obj1.(name)(i,:)), mean(obj2.(name)(i,:)), ...
            std(obj1.(name)(i,:)), std(obj2.(name)(i,:)));
          y.(name).prc68(i,:) = [y.(name).ave(i)-y.(name).sd(i), y.(name).ave(i)+y.(name).sd(i)];
        end
      end
    end

%%
    function y = wtMean(varargin)
      if mod(nargin,2) ~= 0, error('Need to input proper number of argument!'); end
      sum_tmp  = zeros(size(varargin{1}));
      sum_wt = zeros(size(varargin{1}));
      for i = 1:2:nargin
        data = varargin{i};  
        wt   = varargin{i+1};
        sum_tmp  = data.*wt + sum_tmp;
        sum_wt   = sum_wt + wt;
      end
      y = sum_tmp./sum_wt;  
    end

%% 
    function [aveWT,sdWT] = wtMean_SD(ave1,ave2,sd1,sd2)
      aveWT = (ave1 * (1/sd1^2) + ave2 * (1/sd2^2))/((1/sd1^2) + (1/sd2^2));
      sdWT  = 1/sqrt(((1/sd1^2) + (1/sd2^2)));
    end

%%
    function y = convertFLossToLife(model)
      y = Loss.convertFLossToLife(model.para.tspan,model.para.floss,model); 
    end

%%
    function y = convertFLossToOHAnomaly(model)
      y = Loss.convertFLossToOHAnomaly(model.para.tspan,model.para.floss,model); 
    end

%%
    function y = convertFLossToLoss(model)
      y = 1./Loss.convertFLossToLife(model.para.tspan,model.para.floss,model);
    end

%%
    function y = convertFLossToOH(model)
      oh_anomaly = Loss.convertFLossToOHAnomaly(model.para.tspan,model.para.floss,model); 
      y = 1/Loss.lifeP12.oh .* (1+oh_anomaly./100);
    end

%% 
    function plotObj = AddCalc(plotObj,nameList,model)
      for i = 1:length(nameList)
        if isempty(plotObj.(nameList(i)))
          if isempty(model.para)
            model.para = plotObj; 
          end
          plotObj.(nameList(i)) = MyCalc.getDeltaBIOorFF(nameList(i),model); 
        end
      end
      clearNameList = ClearMemory.getFieldList(plotObj,nameList);
      for name = clearNameList, plotObj.(name) = []; end
      if isempty(plotObj.life) && ~isempty(plotObj.floss)
        plotObj.life = Loss.convertFLossToLife(plotObj.tspan,plotObj.floss,model); 
      end
    end

%% 
    function y = getDeltaBIOorFF(name,model)
      switch name
        case 'd13Cbio'
          ems = Emission().getCH4only(model,["anth_bio","natr_bio"]);
          y = MyCalc.wtMean(model.para.d13Canth_bio, ems.anth_bio, model.para.d13Cnatr_bio, ems.natr_bio);
        case 'd13Cff'
          ems = Emission().getCH4only(model,["anth_ff","geo"]); 
          y = MyCalc.wtMean(model.para.d13Canth_ff, ems.anth_ff, model.para.d13Cgeo, ems.geo);
        case 'dDbio'
          ems = Emission().getCH4only(model,["anth_bio","natr_bio"]); 
          y = MyCalc.wtMean(model.para.dDanth_bio, ems.anth_bio, model.para.dDnatr_bio,ems.natr_bio);
        case 'dDff'
          ems = Emission().getCH4only(model,["anth_ff","geo"]); 
          y = MyCalc.wtMean(model.para.dDanth_ff, ems.anth_ff, model.para.dDgeo, ems.geo);
        case 'life' 
          y = MyCalc.convertFLossToLife(model);
        otherwise
          error('No such a fieldname in para!')
      end
    end

%% 
    function classType = ConvertPlotToClassType(plotType)
      switch plotType
        case {'Atmosphere'}
          classType = "Atmosphere";
        case {'ECH4'}
          classType = "Emission";
        case {'Parameter','Isotope'}
          classType = "ParameterSourceSink";
        case {'SourceFraction'}
          classType = "SourceFraction";
        otherwise
          classType = "Struct"; 
      end  
    end     

%% 
    function plotObj = convertCalcType(calcType,plotObj,nameList,periodList) 
      if isa(plotObj.(nameList(1)),'struct'), return; end  
      if isa(periodList,'cell') && length(periodList) == 1 
        periodList = periodList{1};
      end

      clearNameList = ClearMemory.getFieldList(plotObj,nameList);
      plotObj = ClearMemory.Select(plotObj,clearNameList);

      for name = nameList 
        switch calcType
          case 'Timeseries'
            return; 
          case 'PeriodMean'
            plotObj.(name) = MyCalc.PeriodMean(plotObj.tspan,plotObj.(name),periodList);
          case 'PeriodDiff'
            if ismember(name,["floss","totLoss","oh"])
              plotObj.(name) = MyCalc.PeriodDiff(plotObj.tspan,plotObj.(name),periodList,"percent"); 
            else
              plotObj.(name) = MyCalc.PeriodDiff(plotObj.tspan,plotObj.(name),periodList);
            end
          case 'LinearTrend'
            plotObj.(name) = MyCalc.LinearTrend(plotObj.tspan,plotObj.(name),periodList,"no");
          case 'LinearTrendDiff' 
            plotObj.(name) = MyCalc.LinearTrendDiff(plotObj.tspan,plotObj.(name),periodList,"no");
          otherwise
            error('No such a type name! Please check the script.')
        end
      end

      if isa(periodList,'cell') 
        plotObj.tspan = mean(mean(reshape(horzcat(periodList{:}),2,2,[]),3),2); 
      else
        plotObj.tspan = mean(periodList,'all');
      end
    end

%%
    function y = PeriodMean(objTspan,obj,periodList)
      if length(objTspan) == 1
        y = obj; 
        return; 
      end 

      y = zeros(size(periodList,1),size(obj,2));

      for i = 1:size(periodList,1)
        if periodList(i,1) == periodList(i,2) 
          indexTspan = objTspan == periodList(i,1);
        else
          indexTspan = floor(objTspan) >= periodList(i,1) & floor(objTspan) <= periodList(i,2);
        end
        y(i,:) = mean(obj(indexTspan,:),1,'omitmissing'); 
      end
    end

%% e.g) y = MyCalc.LinearTrend(paraBest.tspan,paraBest.life,[2000 2015],'yes')
    function [y, slope, fit] = LinearTrend(objTspan,obj,tspanTarget,isPlotFig)
      p = zeros(2,size(obj,2));
      
      if isnan(tspanTarget)
        indexTspan = true(size(objTspan));
      else
        indexTspan = floor(objTspan) >= tspanTarget(1) & floor(objTspan) <= tspanTarget(2);
      end

      f = zeros(size(obj(indexTspan,:)));
      fit.tspan = objTspan(indexTspan).';

      for i = 1:size(obj,2)
        p(:,i) = polyfit(fit.tspan,obj(indexTspan,i),1); 
        f(:,i) = polyval(p(:,i),fit.tspan);
      end

      y = p(1,:); %linear trend against ensembles

      slope.ave = mean(y,2);
      slope.med = prctile(y,50,2);
      slope.p68 = prctile(y,[2.5 97.5],2);
      slope.p95 = prctile(y,[2.5 97.5],2);

      fit.ave = mean(f,2);
      fit.p95 = prctile(f,[2.5 97.5],2);

      if strcmp(isPlotFig,"yes")
        figure; plot(fit.tspan,obj(indexTspan,:),'-','Color',[0.7,0.7,0.7],'LineWidth',0.5); hold on
        plot(fit.tspan,fit.ave,'-','Color','r','LineWidth',2);
        plot(fit.tspan,fit.p95,'--','Color','r','LineWidth',1.5);
      end
    end

%%
    function [y, slope, fit] = LinearTrendDiff(objTspan,obj,tspanTarget,isPlotFig)
      p = zeros(2,size(obj,2),2);
      
      for j = 1:2
        if size(tspanTarget,2) ~=2
          error('Need to set size(tspanTarget,2) ==2 for LinearTrendDiff')
        else
          indexTspan = floor(objTspan) >= tspanTarget(1,j) & floor(objTspan) <= tspanTarget(2,j);
        end
        
        f = zeros(size(obj(indexTspan,:)));
        fit.tspan = objTspan(indexTspan).';
        
        for i = 1:size(obj,2)
          p(:,i,j) = polyfit(fit.tspan,obj(indexTspan,i),1);
          f(:,i,j) = polyval(p(:,i),fit.tspan);
        end
      end

      y = p(1,:,2) - p(1,:,1);  %linear trend against ensembles
      slope.ave = mean(y,2);
      slope.med = prctile(y,50,2);
      slope.p68 = prctile(y,[2.5 97.5],2);
      slope.p95 = prctile(y,[2.5 97.5],2);

      fit.ave = mean(f,2);
      fit.p95 = prctile(f,[2.5 97.5],2);

      if strcmp(isPlotFig,"yes")
        figure; plot(fit.tspan,obj(indexTspan,:),'-','Color',[0.7,0.7,0.7],'LineWidth',0.5); hold on
        plot(fit.tspan,fit.ave,'-','Color','r','LineWidth',2);
        plot(fit.tspan,fit.p95,'--','Color','r','LineWidth',1.5);
      end
    end

%%
function y = PeriodDiff(objTspan,obj,periodList,varargin)
      if isa(periodList,'cell') 
        peropdListHorz = horzcat(periodList{:});
        tspanTarget1 = reshape(peropdListHorz(1,:).',2,[]).'; %1st period
        tspanTarget2 = reshape(peropdListHorz(2,:).',2,[]).'; %2nd period
      else
        tspanTarget1 = periodList(1,:);
        tspanTarget2 = periodList(2,:);
      end

      objTarget1 = MyCalc.PeriodMean(objTspan,obj,tspanTarget1);
      objTarget2 = MyCalc.PeriodMean(objTspan,obj,tspanTarget2);

      y = objTarget2 - objTarget1;

      if ~isempty(varargin) 
        switch varargin{1}
          case "percent", y = y./objTarget1.*100; %oh or totLoss
        end
      end
    end

 %% 
    function y = infoType(type,period)
      switch type
        case 'LinearTrend'
          y = strcat(type," for ",num2str(period(1)),"-",num2str(period(2)));
        case 'PeriodDiff'
          y = strcat(type," ",num2str(period(2,1)),"-",num2str(period(2,2)), ...
            " minus ",num2str(period(1,1)),"-",num2str(period(1,2)));
        case 'PeriodMean'
          y = strcat(type," for ",num2str(period(1)),"-",num2str(period(2)));
      end
    end  

%% 
    function y = Anomaly(AnomSet,tspan,y) 
      targetMean = mean(interp1(tspan,y,AnomSet.period(1):0.1:AnomSet.period(2)),1);
      assert(~any(isnan(targetMean)),"targetMean contains NaN! Please check the input data.")

      switch AnomSet.type
        case 'absolute'
          y = y - targetMean;
        case 'percent'
          y = (y - targetMean)/targetMean * 100;
      end  
    end

  end
end


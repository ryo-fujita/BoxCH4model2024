% Adjust time axis for timeseries plot between obj.Lim(1) and obj.Lim(2)
% For example, see Fig. 1 and Fig. 2 in Fujita et al. JGR
% Ryo Fujita 2024
classdef TspanScale
  properties 
    flag = 'on'; %'on' or 'off'
    tspan_org
    plotPeriod 
  end

  properties(Dependent)
    Lim  
    Grid 
    tspan
    scale
  end

  methods
    function obj = TspanScale(varargin)
      if nargin > 0
        obj.tspan_org  = varargin{1};
        obj.plotPeriod = varargin{2};
      end
    end

%%
    function y = get.Lim(obj)
      if obj.plotPeriod(1) <= 1750
        y = [1700 1950];  
      elseif obj.plotPeriod(1) == 1800
        y = [1800 1950];
      else
        y = [1850 1950];
      end
    end

%%
    function y = get.Grid(obj)
      if obj.Lim(1) == 1700 && obj.Lim(2) == 1950
        y = 50;  
      elseif obj.Lim(1) == 1800 && obj.Lim(2) == 1950 
        y = 25;
      elseif obj.Lim(1) == 1850 && obj.Lim(2) == 1950 
        y = 25;
      end
    end

%%
    function y = get.tspan(obj)
      if isempty(obj.tspan_org), error('Need to set tspan_org!'); end
        if strcmp(obj.flag,'on') 
          y = obj.convertTspanScaleDown;
        else
          y = obj.tspan_org;
        end
    end

%%
     function y = get.scale(obj)
       y = obj.Grid/(obj.Lim(2)-obj.Lim(1));
     end

%%
    function SetXlimSubplot(obj,xLim)
      if strcmp(obj.flag,'on') && xLim(1) <  obj.Lim(2)
        obj.SetTspanAxisScaled(xLim);
        xlim([obj.convertTspanScaleDown]);
      else
        xlim(xLim);
        xticks('auto');
      end
    end

%%
    function tspanNew = convertTspanScaleDown(obj)
      tspanNew    = zeros(size(obj.tspan_org));
      tspanNew_st = obj.Lim(2) -  obj.Grid;
      
      for i = 1:length(obj.tspan_org)
        t_org = obj.tspan_org(i);
        if t_org <= obj.Lim(1) 
          tspanNew(i) = tspanNew_st - (obj.Lim(1)-t_org);
        elseif t_org <= obj.Lim(2)
          tspanNew(i) = tspanNew_st + (t_org - obj.Lim(1))*obj.scale;
        else
          tspanNew(i) = t_org;
        end
      end
    end

%%
    function tspanOrg = invertTspanOrg(obj,tspanNew)
      tspanNew_st = obj.Lim(2) -  obj.Grid;
      tspanOrg    = zeros(size(tspanNew));
      for i = 1:length(tspanNew)
        if tspanNew(i) <= obj.Lim(2)
          tspanOrg(i) = obj.Lim(1) + (tspanNew(i) - tspanNew_st)/obj.scale;
        else
          tspanOrg(i) = tspanNew(i);
        end
      end
    end

%% 
    function SetTspanAxisScaled(obj,xLimOrg)
      if obj.Lim(1) == 1700 && obj.Lim(2) == 1950 
        XtickNew       = [1900, 1910, 1920, 1930, 1940, 1950, 2000];
        XTickLabelNew  = {'','1750','','1850','','1950','2000'};

      else
        XtickNew = obj.tspanTickNew(xLimOrg);
        XTickLabelNew = obj.tspanTickNewLabel(xLimOrg);
      end

      set(gca,'box','on','Xtick',XtickNew,'XTickLabel',XTickLabelNew,'FontSize',12)
    end

%%
    function y = tspanTickNew(obj,xlim_org)
      tspanNew_st = obj.Lim(2) -  obj.Grid;
      if obj.Lim(1) == 1700 
       y = [tspanNew_st:obj.Grid^2/diff(obj.Lim):tspanNew_st+obj.Grid-obj.Grid^2/diff(obj.Lim) ...
            obj.Lim(2):obj.Grid/5:xlim_org(2)];
      elseif obj.Lim(1) == 1850 || obj.Lim(1) == 1800
       y = [tspanNew_st:obj.Grid^2/diff(obj.Lim):tspanNew_st+obj.Grid-obj.Grid^2/diff(obj.Lim) ...
            obj.Lim(2):obj.Grid:xlim_org(2)]; 
      end
    end

%%
    function y = tspanTickOrg(obj,tspanTickNew)
      for i = 1:length(tspanTickNew)
        y(i) = xlim_org(1):obj.Grid:xlim_org(2);
      end
    end

%%
    function y = tspanTickNewLabel(obj,xlim_org)
      ax = gca;
      tspanTickNew = obj.tspanTickNew(xlim_org);
      tspanTickOrg = obj.invertTspanOrg(tspanTickNew);
      y = cell(1,length(tspanTickNew));

      for i = 1:length(tspanTickNew)
        if i == 1 && xlim_org(1) == obj.Lim(1)
          y(i) = cellstr(num2str(xlim_org(1)));
        elseif i > 1 && xlim_org(1) ~= obj.Lim(1) && xlim_org(1) <= tspanTickOrg(i) ...
        && mod(tspanTickOrg(i),obj.Grid) == 0 && isempty(y{i-1})
          y(i) = cellstr(num2str(tspanTickOrg(i)));
        elseif tspanTickNew(i) < obj.Lim(2)
          y(i) = cellstr('');
        elseif mod(tspanTickOrg(i),obj.Grid) ~= 0
          yyyy = cellstr(num2str(tspanTickNew(i)));
          y(i) = extractBetween(yyyy{1},3,4);
        else
          y(i) = cellstr(num2str(tspanTickNew(i)));
        end
      end
    end

%%
    function y = tspanTickNewLabel_org(obj,xlim_org)
      tspanTickNew = obj.tspanTickNew(xlim_org);
      y = cell(1,length(tspanTickNew));
      for i = 1:length(tspanTickNew)
        if i == 1 
          y(i) = cellstr(num2str(xlim_org(1)));
        elseif tspanTickNew(i) < obj.Lim(2)
          y(i) = cellstr('');
        else
          y(i) = cellstr(num2str(tspanTickNew(i)));
        end
      end
    end
  end

end % end class
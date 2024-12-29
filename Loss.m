% Handles loss rates in total CH4, 13CH4, CH3D, and 14CH4.
% Ryo Fujita 2024

classdef Loss
  properties
    tspan 
    CH4  
    C13  
    CHD  
    C14  
  end
 
  properties(Constant)
    lifeP12 = struct('tot',9.1,'oh',11.18,'soil',150,'strat',120,'cl',200); % Prather et al. (2012), GRL
  end
 
  methods
%%
    function obj = Loss(model) 
      if nargin == 0, return; end 
 
      obj(1,size(model.para,2)) = obj;

      for i = 1:size(model.para,2)
        obj(i).tspan = model.para(i).tspan;
      end

      obj(1) = SetLoss(obj(1),model,model.para(1));
    end
 
%%
    function obj = SetLoss(obj,model,para) 
      obj.CH4  = obj.SetLossCH4("CH4",model,para);
      obj.C13  = obj.SetLossCH4("C13",model,para);
      obj.CHD  = obj.SetLossCH4("CHD",model,para);
      obj.C14  = obj.SetLossCH4("C14",model,para);
    end
 
%%
    function out  = SetLossCH4(obj,gasID,model,para)
      if isempty(para.life) && ~isempty(para.floss)
        para.life = Loss.convertFLossToLife(obj.tspan,para.floss,model);
      end
 
      switch gasID
        case 'CH4'
          out(1).tot = 1./para.life;
        case 'C13'
          out(1).tot = obj.CH4.tot .* 1./para.KIEC; 
        case 'CHD'
          out(1).tot = obj.CH4.tot .* 1./para.KIED; 
        case 'C14'
          loss_R      = MassConstants.loss_R;
          out(1).tot = obj.C13.tot .* 1./para.KIEC + loss_R; 
      end
    end  
 
%%
    function y = isFinite(obj)
      y = all(isfinite(obj.CH4.tot),1) & all(isfinite(obj.C13.tot),1) & ...
        all(isfinite(obj.CHD.tot),1) & all(isfinite(obj.C14.tot),1);
    end
 
  end
 
%%
  methods (Static) 
    function y = convertFLossToLife(tspan,floss,model)
      loss_tot = 1./Loss.getLifeOH(tspan,model) + 1./Loss.lifeP12.soil + ...
        1./Loss.lifeP12.strat + 1./Loss.lifeP12.cl;  
      y = 1./(floss.*loss_tot);
    end

%%
    function [yr_oh,ohAnomaly] = readOHAnomaly(model)
      switch model.lossOH
        case "Stevenson20"
          [yr_oh,ohAnomaly] = Loss.readOHAnomalyStevenson2020;
        case "P12_const"
          yr_oh = [0;2200]; 
          ohAnomaly = zeros(2,1); 
        otherwise
          error(strcat("No such  lossOH exists!: ",model.lossOH)) 
      end
    end

%%
    function y = getLifeOH(tspan,model)
      [yr_loss,ohAnomaly] = Loss.readOHAnomaly(model);

      if size(tspan,1) == 1
        tspan = tspan.'; 
      end

      y = Loss.lifeP12.oh./(1+interp1(yr_loss,ohAnomaly,tspan)./100);

      if tspan(end) > yr_loss(end) && any(isnan(y))
        life_end = Loss.lifeP12.oh./(1+ohAnomaly(end)./100);
        y(isnan(y)) = life_end.* ones(sum(isnan(y)),1);
      end
    end

%%
    function [yr_oh,ohAnomaly] = readOHAnomalyStevenson2020()
      [yr_oh,ohAnomaly] = MyFunc.readTextFreefmt("./loss/ohAnomaly_Stevenson2020.txt", '%n %n', '\t');
    end

%%
    function y = convertFLossToOHAnomaly(tspan,floss,model)
      [yr_loss,ohAnomaly] = Loss.readOHAnomaly(model);
      if size(tspan,1) == 1, tspan = tspan.'; end
      y = (floss.*(1+interp1(yr_loss,ohAnomaly,tspan)./100)-1).*100;
    end

%% 
    function y = shortName()
      y = "loss";
    end
  end
end
 

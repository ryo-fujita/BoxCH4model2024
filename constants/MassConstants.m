% A class that aggregates the constants properties. 
% Located in constants/.
% Ryo Fujita 2024

classdef MassConstants
  properties (Constant)
    % VPDB scale (Lassey et al., 2007)
    R13C_std = 0.0112372 

    % VSMOW scale https://en.wikipedia.org/wiki/Vienna_Standard_Mean_Ocean_Water
    RD_std   = 155.76*10^(-6) 

    % Bq gC-1 (Lassey et al., 2007; Stuiver, 1980) 
    AabsPergC= 0.2260 

    % mole(14C)
    Bq       = 433.2 * 10^(-15) 

    % Tg(CH4)/ppb Prather et al. 2012 2021.3.16 apply
    factor   = 2.75 

    %  gram per mole for total CH4
    W_CH4    = 16.01

    % radiocarbon decay constant (yr-1) 
    loss_R   = 1/8267 
  end

  properties (Dependent)
    Aabs
    W_C13
  end
  
  methods
    function a = get.Aabs(obj) %1gC -> 1/12 mole(total C)
      a  = obj.AabsPergC * obj.Bq /(1/12.01);
    end

    function a = get.W_C13(obj)
      a  = 17/16 * obj.R13C_std * (1 -47/1000);
    end

  end
end


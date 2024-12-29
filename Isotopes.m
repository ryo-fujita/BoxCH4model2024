% A class that aggregates the functions for handling isotope signatures 
% and related emission calculations
% Ryo Fujita 2024

classdef Isotopes 
  properties
    d13C
    dD
    D14C
    P14C
    d13C_org
    dD_org
  end

  properties(Constant)
    emsPWRFile = "./emission/PWRdata_Lassey1960-1971_IAEA1972-2016.txt"
    D14CH4bioFileNameHist = "./emission/D14CH4bio_historical_1750-2015.txt"
  end
  
  methods(Static)
    function y = cnameUse(cname,model,para)
      paraName = model.nameList.para;

      switch cname
        case 'bb'
          y = cname; 
        case {'rumi','rice','wast','wet','trmt','anim'}
          if (sum(ismember(paraName,'fbio'),'all') == 1 ...
            && ~isempty(para.fbio)) ...
            || sum(ismember(paraName,'d13Cbio'),'all') == 1 
            y = 'bio';
          else
            switch cname
              case {'rumi','rice','wast'}
                if sum(ismember(paraName,'fanth_bio'),'all') == 1 ...
                  && ~isempty(para.fanth_bio)
                  y = 'anth_bio';
                else
                  error('Wrong parameter name or isotope setting for anth_bio!')
                end
              case {'wet','trmt','anim'}
                if sum(ismember(paraName,'fwet'),'all') == 1 ...
                  && ~isempty(para.fwet) && strcmp(cname,'wet')
                  y = 'wet';
                elseif sum(ismember(paraName,'fnatr_bio'),'all') == 1 ...
                  && ~isempty(para.fnatr_bio)
                  y = 'natr_bio';
                else
                  error('Wrong parameter name or isotope setting for natr_bio!')
                end
            end
          end
        case {'gas','coal','rco','otherff'}
          if (sum(ismember(paraName,'fff'),'all') == 1 && ...
            ~isempty(para.fbio)) || sum(ismember(paraName,'d13Cff'),'all') == 1
            y = 'ff';
          elseif sum(ismember(paraName,'d13Canth_ff'),'all') == 1
            y = 'anth_ff';
          end
        case 'geo'
          if (sum(ismember(paraName,'fff'),'all') == 1 && ...
            ~isempty(para.fbio)) || sum(ismember(paraName,'d13Cff'),'all') == 1
            y = 'ff';
          elseif sum(ismember(paraName,'d13Cgeo'),'all') == 1
            y = 'geo';
          end
        otherwise
          y = cname;
      end
    end

%% 
    function [E_CXX,delta] = ReadIsotopeBioWithTau(...
      gasID,E_CH4,d13Cs,model,para,cname) 
      if isempty(d13Cs)
        error(strcat("d13C",cname," is empty!")); 
      end

      numData  = model.numData; 
      numCase  = model.numCase; 

      paraName = model.nameParaList(1).all.';
      tspanEms = para.tspan.';  

      M        = MassConstants;
      Aabs     = M.Aabs; 

      if sum(ismember(paraName,'tau'),'all') == 1 
        tau = para.tau; 
      else
        switch cname
        case {'rumi','rice','wast'}
          tau     = 1 * ones(numData,numCase);
        case {'wet','trmt','anim'}
          tau     = 20 * ones(numData,numCase);
        case 'bb'
          tau     = 6 * ones(numData,numCase);
        end
      end

      tauList = [0.1 0.5 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

      if ~all(tau >= 0.1 & tau <= 20,'all')
        disp('para.tau contains out of ranges'); 
      end

      deltaBIO = Isotopes.SetDeltaBIO(model,tspanEms); 
      delta  = zeros(numData ,numCase);
      for i=1:length(tspanEms)
        delta(i,:) = interp1(tauList,deltaBIO(i,:),tau(i,:));
      end
      E_CXX  = E_CH4 .* Aabs.*(1 + delta./1000)./(0.975./(1+d13Cs./1000)).^2;
   end

%%
    function y = EpsilonBIO(obj,gasID,cname,deltaBIO)
      switch gasID
        case 'C13'
          y = obj.d13C_org.(cname).def - deltaBIO;
      end
    end

%% 
    function DeltaBIO = SetDeltaBIO(model,tspanPara)
      [tspanBIO,DeltaBIOWith22tau] = Isotopes.SetDeltaBIOtauList(model); 

      t1 = tspanPara < model.yrInitial-0.5;
      t2 = tspanPara >= model.yrInitial-0.5 & tspanPara <  tspanBIO(end);
      t3 = tspanPara >= tspanBIO(end);
      i1 = find(t1); 
      i2 = find(t2); 
      i3 = find(t3); 

      if(model.yrInitial-0.5 < tspanBIO(1))
        disp("DeltaBIO is extrapolated!"); 
      end

      if sum(i1) >= 1 
        DeltaBIO(i1,:) = interp1(...
          tspanBIO, DeltaBIOWith22tau, model.yrInitial-0.5, 'linear').*ones(length(i1),1); 
      end 

      if sum(i2) >= 1
        DeltaBIO(i2,:) = interp1(tspanBIO, DeltaBIOWith22tau, tspanPara(i2), 'linear'); 
      end
     
      if sum(i3) >= 1
        DeltaBIO(i3,:) = DeltaBIOWith22tau(end,:).*ones(length(i3),1); 
      end
    end 

%% 
    function [tspan,DeltaBIOWith22tau] = SetDeltaBIOtauList(model)
      D14CH4bioFiles = Isotopes.getD14CH4bioFileNameRead(model); 

      for filename = D14CH4bioFiles
        if exist(filename,'file') ~= 0 
          if(strcmp(filename,D14CH4bioFiles(1)))
            y = MyFunc.readText(filename,23,false);
          else
            y = unique([y; MyFunc.readText(filename,23,false)]);
          end
        else
          error(strcat(filename, ...
            " does not exist. You need to prepare it before main simulation."))
        end
      end

      tspan = y(:,1);
      DeltaBIOWith22tau = y(:,2:23);
    end

%%
    function y = getD14CH4bioFileNameRead(model)
      switch(model.issp)
        case 0
          if model.yrFinal <= 2015
            y = Isotopes.D14CH4bioFileNameHist;
          else
            error(strcat("There is no applicable D14Cbio file for issp = 0 after ", ...
              num2str(model.yrFinal)))
          end
        otherwise
          if model.yrInitial < 2015 
            y = [Isotopes.D14CH4bioFileNameHist, Isotopes.D14CH4bioFileNameSSP(model.issp)];
          elseif model.yrInitial >= 2015
            y = Isotopes.D14CH4bioFileNameSSP(model.issp);
          end
      end
    end

%% 
    function y = getD14CH4bioFileNameWrite(model)
      if (model.yrInitial == 1750 && model.yrFinal == 2015 && model.issp == 0)
        y = Isotopes.D14CH4bioFileNameHist;
      elseif(model.yrInitial == 2015 && model.yrFinal == 2100 && model.issp ~= 0)
        y = Isotopes.D14CH4bioFileNameSSP(model.issp);
      else
        error("The model.yrInitial and/or model.yrFinal are " + ...
          "not applicable for D14Cbio calculation!") 
      end
    end

%%

    function y = D14CH4bioFileNameSSP(issp)
      y = strcat("./emission/D14CH4bio_",MyNameList.ssp(issp),"_2015-2100.txt"); 
    end

%% 
    function E_C14 = ComputeEmissionNPR(model,para) 
      numData  = model.numData; 
      numCase  = model.numCase;  
      E_C14    = zeros(numData,numCase);
      M        = MassConstants;
      Bq       = M.Bq;
      W_CH4    = M.W_CH4;
      
      if size(para.tspan,1) == 1
        tspanEms = para.tspan.';  
      else
        tspanEms = para.tspan; 
      end 

      [tspan_PWR,PWRdata] = Isotopes.SetPWRInventorydata;
      phi = para.phi * Bq * 10^9 ;% GBq(14CH4)*GWe-yr^-1 -> mole(14C)*GWe-yr^-1 : 

      t1 = tspanEms < 1960.5;
      t2 = 1960.5 <= tspanEms & tspanEms <= tspan_PWR(end); 
      t3 = tspan_PWR(end) < para.tspan; 
      i1 = find(t1);
      i2 = find(t2);
      i3 = find(t3); 

      yday = ones(numData,1).*365 + (mod(fix(tspanEms),4) == 0); % vector of days in each year

      if sum(i1) >= 1
        E_C14(i1,:) = zeros(length(i1),1) .*  phi(t1,:) ; 
      end
 
      if sum(i2) >= 1
        E_C14(i2,:) = (interp1(tspan_PWR,PWRdata,tspanEms(t2),'linear')./ ...
          (24 .* yday(t2)) .* W_CH4 .* 10^(-12)) .* phi(t2,:); 
      end  

      if sum(i3) >= 1
        E_C14(i3,:) = PWRdata(end)./(24 .* yday(t3)) .* W_CH4 .* 10^(-12) ...
          .* phi(t3,:); 
      end
    end

%% 
    function [tspan_PWR,PWRdata] = SetPWRInventorydata()
      for filename = Isotopes.emsPWRFile
        if exist(filename,'file') ~= 0 
          switch(extractAfter(filename,"emission/"))
            case {"PWRdata_Lassey1960-1971_IAEA1972-2016.txt"}
              [tspan_PWR_tmp,PWRdata_tmp] = MyFunc.readTextFreefmt(filename, '%n %n', '\t');
            otherwise
              error(strcat(filename, " is not used for the current simulation."))
          end

          if(strcmp(filename,Isotopes.emsPWRFile(1)))
            tspan_PWR = tspan_PWR_tmp;
            PWRdata   = PWRdata_tmp;            
          else
            tspan_PWR = [tspan_PWR; tspan_PWR_tmp];
            PWRdata   = [PWRdata; PWRdata_tmp];    
          end
        else
          error(strcat(filename," does not exist. ..." + ...
            "You need to prepare it before main simulation.")) 
        end        
      end
    end

%% 
    function [delta, E_CXX] = computeDeltaBIOWithTau_Integral(...
      obj,gasID, cname, E_CH4, iData, yr, tau, model)  % all input variable are set to be scalar
      issp     = model.issp;
      M        = MassConstants;
      R13C_std = M.R13C_std; 
      Aabs     = M.Aabs; 
      
      file.name1 = ...
        './emission/d14CO2&d13CO2_Graven2017_Reimer2013add_20220112.txt';
      file.name2 = ...
        './emission/d14CO2&d13CO2_ssp_Graven2018.txt';

      [inpdata.yr_CO2, inpdata.d14CO2, inpdata.d13CO2] ...
        = MyFunc.readTextFreefmt(file.name1, '%n %n %n', '\t');   
                         
      A = importdata(file.name2, '\t', 3);
      inpdata.yr_CO2_ssp = A.data(:,1); %SSP119  SSP126	SSP245	SSP3B	SSP534os	SSP5B

      [inpdata.d14CO2_ssp] = ...
        horzcat(A.data(:,2),A.data(:,4),A.data(:,6),A.data(:,8),...
        A.data(:,10),A.data(:,12));
      [inpdata.d13CO2_ssp] = horzcat(A.data(:,3),A.data(:,5),A.data(:,7),...
        A.data(:,9),A.data(:,11),A.data(:,13));

      switch gasID
        case 'C13'
          loss       = 0;
          yr_atm     = inpdata.yr_CO2;
          d_atm      = inpdata.d13CO2;
          if yr > 2015.5 && issp ~= 0
            yr_atm_ssp = inpdata.yr_CO2_ssp;
            d_atm_ssp  = inpdata.d13CO2_ssp(:,issp);
          end
        case 'C14'
          loss       = MassConstants.loss_R;
          d13Cs      = obj.d13C.(cname)(iData,:);
          yr_atm     = inpdata.yr_CO2;
          d_atm      = inpdata.d14CO2;
          if yr > 2015.5 && issp ~= 0
            yr_atm_ssp = inpdata.yr_CO2_ssp;
            d_atm_ssp  = inpdata.d14CO2_ssp(:,issp);
          end
      end 

      xmax  = 200; %maximum integral year
      fun = @(tt) (1 + interp1(yr_atm,d_atm,yr - tt,'linear')./1000) ...
        .*exp(-tt./tau)./tau.*exp(-loss.*tt); 
      
      if yr - xmax < d_atm(1)
        error(strcat('You need to prepare d14CO2 data before', num2str(yr - xmax))); 
      end

      if yr <= 2015.5
        P_src = integral(fun,0,xmax,'ArrayValued',true); 
      else
        fun_ssp  = @(tt) (1 + interp1(yr_atm_ssp,d_atm_ssp, yr - tt,'linear')./1000) ...
          .*exp(-tt./tau)./tau.*exp(-loss.*tt);
        if yr - xmax < 2015.5   
          P_src = integral(fun_ssp,0,yr-2015.5,'ArrayValued',true) + ...
            integral(fun, yr-2015.5, xmax, 'ArrayValued', true);
        else
          P_src = integral(fun_ssp, 0, xmax, 'ArrayValued', true);
        end  
      end

      delta  = (P_src - 1).*1000;

      switch gasID
        case 'C13'
          E_CXX = E_CH4 .* R13C_std .* ...
            (1 + delta./1000) ./(R13C_std.*(1 + delta./1000) + 1);
        case 'C14'
        % E_C14_{S} = E_CH4 * A_{S}  (measured sample activity)
          E_CXX = E_CH4 .* Aabs.*(1 + delta./1000)./(0.975./(1+d13Cs./1000)).^2;
      end  
    end 
  end
end

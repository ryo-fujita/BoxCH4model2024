% Handles the calculations for atmospheric CH4 concentration, 
% d13C, dD, and D14C.
% Ryo Fujita 2024

classdef Atmosphere 
  properties
    termName 
    tspan 
    target 
    nameList 
    concNameList 
    CH4    
    d13C   
    dD     
    D14C   
    C13    
    CHD    
    C14    
  end

  methods
%% 
    function obj = Atmosphere(model,varargin)
      if nargin ~= 0
        if model.isLagSmoother || model.isPosterior || model.isPosteriorOSSE || ...
           model.isReconstruction || model.isRunWithoutFiltering
          obj(1).nameList = model.nameList.atmos; 
          obj(1).concNameList = model.nameList.atmosConc; 
          obj = obj.DoSpecialCase(model,varargin{:}); 
          return 
        end

        obj(1,model.numNode) = obj;
        for i = 1:model.numNode
          obj(i).nameList = model.nameList.atmos; 
          obj(i).concNameList = model.nameList.atmosConc; 
          switch i
            case 1
              obj(i).termName = "initial";
              obj(i).tspan    = model.atmosTargetList(1).tspan(1);
              obj(i).target   = model.atmosTargetList(1);
            case model.numNode
              obj(i).termName = "final";
              obj(i).tspan    = obj(i-1).tspan(end):1:model.yrFinal;
            otherwise
              obj(i).termName = model.nameList.term(i-1);
              obj(i).tspan    = model.atmosTargetList(i-1).tspan(1):1:model.atmosTargetList(i-1).tspan(2); 
              obj(i).target   = model.atmosTargetList(i-1);
          end
        end
        obj(1) = obj(1).SetAtmosphereInitial(model);
      end
    end

%% 
    function obj = DoSpecialCase(obj,model,varargin)
      if model.isLagSmoother 
        obj = Smoothing.loadPosterior("atmos",model);
        obj = MultiplyMassFactorToMoleFraction(obj);

      elseif model.isPosterior || model.isPosteriorOSSE 
        model.lagFilterTerm = model.numNode;
        obj = Smoothing.loadPosterior("atmos",model);

      elseif model.isReconstruction 
        obj = model.loadFilteredData("atmos").atmos; 

      elseif model.isRunWithoutFiltering 
        obj(1).tspan  = model.atmosTargetList(1).tspan(1):1:model.yrFinal;
        obj(1).target = model.atmosTargetList(1);
        obj(1) = obj(1).SetAtmosphereInitial(model);
      end

      if ~isempty(varargin) %e.g. obj = obj.ConvertToIsotopeMoleFraction
        for i = 1:length(varargin)
          funcName = varargin{i}; 
          obj = obj.(funcName);
        end
      end
    end

%%
    function obj = SetAtmosphereInitial(obj,model) 
      if model.isPosteriorMeanOnly 
        atmosPost = Atmosphere(Model(model.workdir,"Posterior",9));
        for gasID = model.nameList.atmos
          obj.(gasID) = interp1(atmosPost.tspan,mean(atmosPost.(gasID),2),model.yrInitial);
        end
        obj = ConvertToIsotopeMoleFraction(obj); 

      elseif model.isPosterior %2024.12.23
        obj = obj.SetAtmosphereFromPostBegin(model);

      elseif ismember(model.currentNode,[1,2]) && model.issp == 0
        obj = obj.SetInitialDistribution(model);

      elseif model.currentNode == 1 && model.issp ~= 0
        if model.currentIterNum == 0
          if model.currentRun==0, iRun=1; else, iRun=model.currentRun; end 
          atmosHist = matfile(strcat(model.workdir,'/objects_SmoothSummary_Full.mat')).atmos;
          atmosHist = atmosHist(iRun,end-1);
        else
          atmosHist = model.atmos(1);
        end
        obj = obj.SetAtmosphereFromPostLast(atmosHist);

      else %set initial from target value
        obj = obj.SetInitialPointFixed(model);
      end
    end

%% 
    function obj = SetAtmosphere(obj,model,input)
      for gasID = ["CH4","C13","CHD","C14","d13C","dD","D14C"]
        obj.(gasID)  = input.(gasID)(end,:) ;
      end
      resIdx = model.judge(model.currentNode-1).resIdx;
      obj = JudgeAtmosToAtmosTarget.updateResults(resIdx,obj);
    end  

%% 
    function obj = SetAtmosphereFromPostBegin(obj,model)
      obj.nameList = model.nameList.atmos;
      obj.concNameList = model.nameList.atmosConc; 
      
      for gasID = ["CH4","d13C","dD","D14C"], obj.(gasID)  = obj.(gasID)(1,:); end
      
      obj = obj.ConvertToIsotopeMoleFraction;
      obj.tspan = obj.tspan(1):1:model.yrFinal;
    end

%% 
    function obj = SetAtmosphereFromPostLast(obj,atmosPost)
      for gasID = ["CH4","d13C","dD","D14C"], obj.(gasID)  = atmosPost.(gasID)(end,:); end
      obj = obj.ConvertToIsotopeMoleFraction;
    end

%% 
    function obj = SetInitialDistribution(obj,model)
      numCase = model.numCase;
      rng(model.rngSetting); %default, rng('shuffle')
      latinArray = lhsdesign(4,numCase);
      CH40  = obj.target.CH4.ave - ...
        obj.target.CH4.sdev*3 + obj.target.CH4.sdev*6 .*latinArray(1,:); %plus minus 3 sigma
      d13C0 = obj.target.d13C.min + (obj.target.d13C.max - obj.target.d13C.min).*latinArray(2,:);
      dD0   = obj.target.dD.min   + (obj.target.dD.max - obj.target.dD.min).*latinArray(3,:);      
      D14C0 = obj.target.D14C.min + (obj.target.D14C.max - obj.target.D14C.min).*latinArray(4,:);
      obj   = obj.SetValue(numCase,CH40,d13C0,dD0,D14C0);
    end

%% 
    function obj = SetInitialPointFixed(obj,model)
      CH40  = mean([obj.target.CH4.min, obj.target.CH4.max]);
      
      if obj.target.d13C.min == -1000 && obj.target.d13C.max == 1000
        d13C0 = mean([-50,-48]);
      else
        d13C0 = mean([obj.target.d13C.min,obj.target.d13C.max]);
      end
      
      if obj.target.dD.min == -1000 && obj.target.dD.max == 1000
        dD0 = mean([-115,-80]);
      else
        dD0 = mean([obj.target.dD.min,obj.target.dD.max]);
      end
      
      if obj.target.D14C.min == -1000 && obj.target.D14C.max == 1000
        D14C0 = mean([-100,100]);
      else
        D14C0 = mean([obj.target.D14C.min,obj.target.D14C.max]);
      end
      
      obj   = obj.SetValue(model.numCase,CH40,d13C0,dD0,D14C0);
    end

%%
    function obj = SetValue(obj,numCase,CH40,d13C0,dD0,D14C0)
      M           = MassConstants;
      factor      = M.factor;
      R13C_std    = M.R13C_std; 
      RD_std      = M.RD_std; 
      Aabs        = M.Aabs; 

      obj.CH4  = SetNumCase(CH40,numCase) .* factor ;

      obj.d13C = SetNumCase(d13C0,numCase);
      R13C0    = R13C_std * (1 + obj.d13C./1000);
      obj.C13  = obj.CH4 .* (R13C0./(R13C0 + 1));

      obj.dD   = SetNumCase(dD0,numCase);
      RD0      = RD_std * (1 + obj.dD/1000);
      obj.CHD  = obj.CH4 .* (RD0./(RD0 + 1));

      obj.D14C = SetNumCase(D14C0,numCase);
      R14C0    = Aabs*(1 + obj.D14C./1000) ./(0.975./(1 + obj.d13C./1000)).^2; 
      obj.C14  = obj.CH4 .* R14C0; 
    end

%%
    function obj = ConvertToIsotopeRatio(obj)
      M           = MassConstants;
      R13C_std    = M.R13C_std; 
      RD_std      = M.RD_std; 
      Aabs        = M.Aabs; 
 
      for name = obj.concNameList 
        if name == "C13" && ~isempty(obj.C13) && all(eq(size(obj.C13),size(obj.CH4)))
          obj.d13C = (obj.C13./(obj.CH4 - obj.C13)./R13C_std - 1) * 1000; 
        end
        
        if name == "CHD" && ~isempty(obj.CHD) && all(eq(size(obj.CHD),size(obj.CH4)))
          obj.dD   = (obj.CHD./(obj.CH4 - obj.CHD)./RD_std - 1) * 1000; 
        end
  
        if name == "C14" && ~isempty(obj.C14)
        if sum(size(obj.CH4) ~= size(obj.d13C) | size(obj.CH4) ~= size(obj.C14) | size(obj.d13C) ~= size(obj.C14)) ~= 0 % for debug
          disp('find error!')
          error([strcat(num2str(size(obj.CH4)),"  ",num2str(size(obj.d13C)),"  ",num2str(size(obj.C14)),"  ",num2str(size(obj.D14C)))])
        else
          obj.D14C = (obj.C14.*(0.975./(1+obj.d13C./1000)).^2./obj.CH4./Aabs -1)*1000; 
        end
        end
      end
    end 

%% 
    function obj = ConvertToIsotopeMoleFraction(obj)
      M           = MassConstants;
      R13C_std    = M.R13C_std; 
      RD_std      = M.RD_std; 
      Aabs        = M.Aabs; 

      for name = obj.nameList
        if name == "d13C" && ~isempty(obj.d13C) 
          R13C0    = R13C_std * (1 + obj.d13C./1000);
          obj.C13  = obj.CH4 .* (R13C0./(R13C0 + 1));
        end
  
        if name == "dD" && ~isempty(obj.dD) 
          RD0      = RD_std * (1 + obj.dD/1000);
          obj.CHD  = obj.CH4 .* (RD0./(RD0 + 1));
        end
  
        if name == "D14C" && ~isempty(obj.D14C)
          R14C0    = Aabs*(1 + obj.D14C./1000) ./(0.975./(1 + obj.d13C./1000)).^2; 
          obj.C14  = obj.CH4 .* R14C0; 
        end
      end
    end 

%% 
    function obj = MultiplyMassFactorToMoleFraction(obj)
      M           = MassConstants;
      factor      = M.factor;
      for gasID = ["CH4","C13","CHD","C14"] 
        if isempty(obj.(gasID)), continue; end
        obj.(gasID) = obj.(gasID) .* factor;
      end
    end 

%% 
    function obj = DevideMassFactorToMoleFraction(obj)
      M           = MassConstants;
      factor      = M.factor;
      for gasID = ["CH4","C13","CHD","C14"]
        if isempty(obj.(gasID)), continue; end
        obj.(gasID) = obj.(gasID) ./ factor;
      end
    end 
  end

  methods(Static)
    function y = shortName()
      y = "atmos";
    end
  end
end

%%
function y = SetNumCase(val,Ncase)
  if length(val) == 1
    y = val .* ones(1,Ncase);
  elseif length(val) == Ncase
    y = val;
  else
    error('not proper vector size of initialValue!')
  end
end




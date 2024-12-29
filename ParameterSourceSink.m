% Handles the calculations for CH4 source and sink parameters
% Ryo Fujita 2024

classdef ParameterSourceSink 
  properties
    tspan  
    termName 
    target 
    fbb
    fbio  
    fanth_bio 
    fnatr_bio 
    fanth_ff
    Egeo
    tau
    phi
    life
    floss 
    totLoss  
    oh       
    KIEC
    KIED
    d13Cbio
    d13Canth_bio 
    d13Cnatr_bio 
    d13Cbb
    d13Cff 
    d13Canth_ff
    d13Cgeo
    dDbio
    dDanth_bio 
    dDnatr_bio 
    dDbb
    dDff  
    dDanth_ff
    dDgeo
    IAVpercent 
    hyper 
  end

  methods
    function obj = ParameterSourceSink(model) 
      if nargin == 0, return; end

      if obj.isSpecialCase(model)
        obj = obj.DoSpecialCase(model); 
        return 
      end

      obj(model.numNode) = obj; 

      if model.useParaTimeVar
        obj = ParaTimeVar(model).SetIAVpercent(obj); 
      end

      if model.useHyperParameter
        obj = HyperParameter(model).setup(obj); 
      end 

      for iNode = 1:model.numNode
        model.currentNode = iNode;
        switch iNode
          case 1
            obj(iNode).termName = "initial";
            obj(iNode).tspan    = model.paraTargetList(1).tspan(1); 
            obj(iNode).target   = model.paraTargetList(1);
            obj(iNode)          = obj(iNode).SetParameter(model); 

          case model.numNode
            obj(iNode).termName = "final";
            obj(iNode).tspan    = obj(iNode-1).tspan(end)-1:1:model.yrFinal+0.5;
            obj(iNode).target   = obj(iNode-1).target;

          otherwise
            obj(iNode).termName = model.nameList.term(iNode-1);
            obj(iNode).tspan    = model.paraTargetList(iNode-1).tspan(1):...
                                  1:model.paraTargetList(iNode-1).tspan(2);
            obj(iNode).target   = model.paraTargetList(iNode-1);
        end
      end  
    end

%% 
    function obj = SetParameter(obj,model,obj_old) 
      switch model.currentNode
      case 1 
        if model.Spup.step == 1 && ismember(model.Spup.type,["noise","Cauchy"])
          obj = obj.computeParameterInitializeNoise(model);
        else
          obj = obj.computeParameterInitialize(model); 
        end

      case 2 
        obj = obj.computeParameterSpinup(model); 

      case model.numNode 
        obj = obj.computeParameterFinalization(model,obj_old); 

      otherwise 
        obj = obj.computeParameterInLoop(model,obj_old);
      end   
    end

%% 
    function obj = computeParameterInitialize(obj,model) 
      numData = model.numData;
      numCase = model.numCase;
      numPara = model.numParaList(1);

      if model.useHyperParameter
        obj.hyper = obj.hyper.SetParameter(model); 
      end
 
      randArray = obj.createRandArray(numCase,numPara.fDynamic,model.randType); 

      iRand = 0;  
      for name = model.nameList.para
        if obj.target.(name).flag ~= 0
          iRand = iRand + 1; 
        end 
        
        if ~model.useDistPrior %default
          if obj.target.(name).flag ~= 0
            obj.(name) = obj.computeParaFullRange(...
              name,randArray(iRand,:),model) .* ones(numData,1); %uniform 
          else             
            obj.(name) = obj.target.(name).def .* ones(numData,numCase); %fixed
          end

        else %optional
          switch obj.target.(name).distPrior
            case "u" %uniform
              obj.(name) = obj.computeParaFullRange(...
                name,randArray(iRand,:),model) .* ones(numData,1); 
            case "f" %fixed
              obj.(name) = obj.target.(name).def .* ones(numData,numCase); 
            case {"n","bn"} %normal, bounded normal
              obj.(name) = obj.getNormrnd(name,model);
              if strcmp(obj.target.(name).distPrior,"bn")
                obj.(name) = obj.replaceBoundaryForBoundedNormalDist(name,model); 
              end
            otherwise
              error(strcat(obj.target.(name).distPrior,": No such a distPrior!"))
          end
        end

      end
    end

%% 
    function obj = computeParameterInitializeNoise(obj,model) 
      paraTarget = load(strcat(model.workdir,'/tmp_paraTarget_Spinup_initial_',...
                        num2str(model.currentRun,'%05.0f'),'.mat')).paraTarget;
      iIter = model.currentIterNum + 1;
      numCaseSpup = length(paraTarget.(model.nameList.para(1)).def);
      paraNoisePercentSp = model.Spup.paraNoisePercent; 

      if numCaseSpup < iIter
        iIter = mod(iIter,numCaseSpup);
        if iIter == 0
          iIter = numCaseSpup; 
        end
      end

      for name = model.nameList.para
        paraTargetLim = [paraTarget.(name).min, paraTarget.(name).max];

        if paraTarget.(name).flag == 0 
          obj.(name) = model.paraTargetList(1).(name).def .* ones(model.numData,model.numCase);
        else
          def_para = paraTarget.(name).def(iIter);
          sd_para = range(paraTargetLim) * paraNoisePercentSp; 
          
          obj.(name) = (obj.distribution(model.Spup.type,model.numCase) ...
            .* sd_para + def_para) .* ones(model.numData,1);
  
          logic = obj.getLogicExceedBoundary(obj.(name),paraTargetLim);

          if any(logic)
            obj.(name)(1,logic) = def_para; 
          end
        end

      end
    end

%% 
    function y = computeParaFullRange(obj,name,randArray,model)
      if model.Spup.exists && model.currentNode == 1 && model.issp == 0 
        %Use min & max from 1st spin-up results
        paraTarget = load(strcat(model.workdir,'/tmp_paraTarget_Spinup_initial_',...
                          num2str(model.currentRun,'%05.0f'),'.mat')).paraTarget; 
      else
        paraTarget = obj.target;
      end

      y = paraTarget.(name).min + (paraTarget.(name).max - paraTarget.(name).min) .* randArray;
    end

%% 
    function obj = computeParameterSpinup(obj,model)
      switch model.status 
        case 2
          for name = model.nameList.para
            obj.(name) = model.para(1).(name) .* ones(model.numData,1); 
          end
        otherwise, error('No such a status is defined during the Spinup!');
      end 
    end

%% 
    function obj = computeParameterInLoop(obj,model,obj_old)
      iNode   = model.currentNode;
      resIdx  = model.judge(iNode-1).resIdx; 
      rng(model.rngSetting); 

      obj.IAVpercent = obj.IAVpercent.SetParameter(model); 
      
      if model.useHyperParameter
        obj.hyper = obj.hyper.SetParameter(model); 
      end

      numPara = model.numParaList(iNode);
      randArray = obj.createRandArray(model.numCase,numPara.fDynamic-numPara.f1,model.randType);
      iRand = 1;

      for name = model.nameList.para
        paraOldRes = obj_old.(name)(end,resIdx);
        switch obj.target.(name).flag 
          case 0 % no time variations
            obj.(name) = paraOldRes.* ones(model.numData,1); 
          case {2,4} % add noise
            obj = obj.predictParaNew(name,paraOldRes,model,randArray(iRand,:)); 
        end 

        if ~ismember(obj.target.(name).flag,[0,1]), iRand = iRand + 1; end 
      end
    end

%% p(t) = p(t-1) + noise
    function obj = predictParaNew(obj,name,paraOld,model,randArray)
      initialRange = obj.target.(name).max - obj.target.(name).min;

      % IAVRange is the range of uniform distribution or variance of Gaussian distribution
      
      % case: small noise (e.g. Egeo)
      if obj.target.(name).flag == 2 
        IAVRange = model.paraIAVpercent.fix/100 * initialRange; 

      % case: normal noise with time variations (default)
      elseif ~isempty(obj.IAVpercent) && model.useParaTimeVar 
        IAVRange = obj.IAVpercent.(name)./100 .* initialRange;

      % case: noise without time variations of maximum IAV range
      elseif ~isempty(obj.IAVpercent) 
        IAVRange = obj.IAVpercent(1,:)./100 .* initialRange;

      % case: prescribed single maximum IAV range
      else
        IAVRange = model.paraIAVpercent.def/100 * initialRange;
      end

      % get noise
      paraNoise = obj.convertIAVRangeToNoise(IAVRange,model,randArray); 
      
      % Add noise with p(t-1)
      paraNew = paraOld + paraNoise; 

      % Replace p(t) when the value exceeds the prescribed boundary values
      paraNew = obj.replaceBoundary(name,paraOld,model,paraNew);

      % Extend timeseries of parameters with p(t)
      obj.(name) = interpParaOldToNew(paraOld,paraNew,obj.tspan);                  
    end

%% Generate noise
    function y = convertIAVRangeToNoise(obj,IAVRange,model,varargin)
      
      if ~isempty(varargin)
        randUniform = varargin{1}; %uniform distribution [0, 1]
      end 

      switch(model.controlPara.newDist)
        case {'Gauss','Cauchy'}
          y = IAVRange .* obj.distribution( ...
          model.controlPara.newDist,model.numCase);

        case 'Uniform'
          y = IAVRange .* obj.distribution('Uniform',randUniform);
      end
    end

%% If the edge of the Gaussian distribution exceeds min,max, replace the value
    function paraNew = replaceBoundary(obj,name,paraOld,model,paraNew) 
      paraTargetLim = [obj.target.(name).min  obj.target.(name).max]; 
      logic = obj.getLogicExceedBoundary(paraNew,paraTargetLim);

      if strcmp(model.controlPara.exBoundary,"paraOld")
        paraNew(:,logic) = paraOld(logic);

      elseif strcmp(model.controlPara.exBoundary,"NaN")
        paraNew(:,logic) = NaN;

      else
        error(model.controlPara.exBoundary + ...
          " is not available setting for replaceBoundary. Please check the code.") 
      end

    end

%% If the edge of the Gaussian distribution exceeds min,max, find the value
%  until values remain within the boundaries remain.
    function paraNew = replaceBoundaryForBoundedNormalDist(obj,name,model)
      paraTargetLim     = [obj.target.(name).min  obj.target.(name).max]; 
      paraNew = obj.(name);
      logic = obj.getLogicExceedBoundary(paraNew,paraTargetLim);

      iWhile = 0;
      while sum(logic) ~= 0 
        iWhile = iWhile + 1; 
        assert(iWhile<20, "logicParaBoundary cannot be converged within 20 trials. ..." + ...
          "Please consider to change the parameter of stdev.")

        model.numCase = sum(logic);
        paraNew(:,logic) = obj.getNormrnd(name,model,logic);

        logic = obj.getLogicExceedBoundary(paraNew,paraTargetLim);
      end
    end

%% 
    function obj = computeParameterFinalization(obj,model,obj_old) 
      resIdx = model.judge(model.currentNode-1).resIdx; 

      numData = model.numData;

      for name = model.nameList.para 
        obj.(name) = obj_old.(name)(end,resIdx).* ones(numData,1);
      end
    end

%%
    function y = isSpecialCase(obj,model)
      y = model.isLagSmoother ...
        || model.isPosterior ...
        || model.isPosteriorOSSE  ...
        || model.isPosteriorMeanOnly ...
        || model.isReconstruction ...
        || model.isRunWithoutFiltering;
    end

%% 
    function obj = DoSpecialCase(obj,model)
      if model.isLagSmoother 
        obj = Smoothing.loadPosterior("para",model); 

      elseif model.isPosterior || model.isPosteriorOSSE
        model.lagFilterTerm = model.numNode;
        obj = Smoothing.loadPosterior("para",model);

      elseif model.isPriorOnly 
        obj = obj.setDefaultInitial(model); 

      elseif model.isPosteriorMeanOnly
        obj = obj.loadPosteriorMean(model); 

      elseif model.isReconstruction 
        obj = model.loadFilteredData("para").para; 

      else
        %Just get empty object
      end

    end

%% 
    function obj = setDefaultInitial(obj,model)
      obj.tspan = model.paraTargetList(1).tspan(1):1:model.yrFinal+0.5;
      for name = model.nameParaList(1).all
        if model.numCase == 1
          obj.(name) = model.paraTargetList(1).(name).def ...
            .* ones(model.numData,model.numCase); 
        else
          obj.target = model.paraTargetList(1);
          randArray = rand(1,model.numCase);
          obj.(name) = obj.computeParaFullRange(name,randArray,model).* ...
            ones(model.numData,1);
        end
      end      
    end

%% 
    function y = getNormrnd(obj,name,model,varargin)
      ave_para = obj.target.(name).def;
      if obj.target.(name).setHyper
        sd_para = obj.hyper.(name);

        if size(sd_para,2) == model.numCase
          y = normrnd(ave_para,sd_para).* ones(model.numData,1);
        else
          logicTrue = varargin{1};
          y = normrnd(ave_para,sd_para(logicTrue)).* ones(model.numData,1);
        end
      else
        sd_para = obj.target.(name).sd;
        y = normrnd(ave_para,sd_para,1,model.numCase).* ones(model.numData,1);
      end
    end

%% 
    function obj = loadPosteriorMean(obj,model)
      obj.tspan = [model.paraTargetList(1).tspan(1):1:model.yrFinal+0.5].'; 
      if exist(strcat(model.workdir,"/para",model.Info.inventory,...
         "_SmoothSummaryPostFull_meanPrc68.csv"),'file')
        tbl = readtable(strcat(model.workdir,"/para",model.Info.inventory,...
              "_SmoothSummaryPostFull_meanPrc68.csv"));
        for name = model.nameList.para
          obj.(name) = interp1(tbl.fyr,tbl.(name+"_mean"),obj.tspan);
        end

      else
        paraPost  = ParameterSourceSink(Model(model.workdir,"Posterior",9));
        for name = model.nameList.para
          obj.(name) = interp1(paraPost.tspan,mean(paraPost.(name),2),obj.tspan);
        end
      end
    end

%% 
    function y = ismemberHyper(obj,name)
      y = any(ismember(fieldnames(obj.target.(name)).',"setHyper")) && obj.target.(name).setHyper;
    end

%% 
    function y = containsIAVpercentEach(obj)
      y = isa(obj.IAVpercent,'ParaTimeVar') && ~isempty(obj.IAVpercent.termName) ...
        && strcmp(extractBefore(obj.IAVpercent.termName,2),"t");
    end

%%
    function y = isWithinPriorRangeIAV(obj,model)
      y = true; 
      for name = model.nameParaList(model.currentNode-1).fDynamic
          y = y & obj.isWithinPriorRange(name);
      end
    end

%% 
    function y = isWithinPriorRange(obj,name)
      paraParticleNew = obj.(name)(size(obj.(name),1),:);
      y = paraParticleNew >= obj.target.(name).min | paraParticleNew <= obj.target.(name).max;
    end

%%
    function y = nameList(obj)
      y = string(fieldnames(obj)).';
      index  = ~ismember(y,["tspan" "termName" "target" "IAVpercent" "hyper"]); 
      y = y(index);      
    end

  end

  methods(Static)
    function y = distribution(dtype,varargin)
      switch dtype
        case {'Gauss','noise'}
          numCase = varargin{1};
          y = randn(1,numCase);
        case 'Uniform'
          randUniform = varargin{1};
          y = -1 + 2*randUniform; %uniform distribution [-1, 1]
        case 'Cauchy' 
          numCase = varargin{1};
          pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',1);
          y  = random(pd,1,numCase);
        otherwise
          error('No such a distribution name!')
      end
    end
%%
    function y = getLogicExceedBoundary(paraNew,paraTargetLim)
      y = paraNew < paraTargetLim(1) | paraNew > paraTargetLim(2);
    end
%%
    function y = computeParameterLassey07(model,name)
      switch name
        case 'phi'
          phi_ave = 286; phi_std = 26;
          y = (phi_ave + phi_std*randn(1,model.numCase)) .* ones(model.numData,1); 
        case 'life'
          life_ave = 8.6; life_std = 0.5;
          y = (life_ave + life_std*randn(1,model.numCase)) .* ones(model.numData,1);
      end
    end

%% 
    function randArray = createRandArray(numCase,numParaRand,randType)
      if numParaRand > 0
        switch randType
        case "latin" %https://jp.mathworks.com/help/stats/lhsdesign.html
          randArray = lhsdesign(numParaRand,numCase); 
        case "rand"
          randArray = rand(numParaRand,numCase); 
        otherwise
          error('model.randType needs to be set properly!')
        end
      else
        randArray = [];
      end       
    end

%%
    function y = shortName()
      y = "para";
    end

  end

end %% End methods (Static)

%% 
function y = interpParaOldToNew(paraOld,paraNew,tspanPara) 
  tspanNew = tspanPara(2:end).';

  % Repeat paraOld at tspanPara(1) & tspanPara(2) in the new prediction
  % this is due to avoid error at the begining of time integration using ode15s
  % for the period of atmos.tspan 
  % (atmos.tspan(1) = tspanPara(1) + 0.5; atmos.tspan(end) = tspanPara(end) - 0.5)
  y = [paraOld; ...
    interp1([tspanNew(1);tspanNew(end)],[paraOld;paraNew],tspanNew)];
end
% The main program called by RunModel.m. It includes basic 
% properties and numerous functions for atmospheric box model. 
% It can also be run independently (e.g., Model().main). 
% Ryo Fujita 2024

classdef Model 
  properties (Constant, Access = public)
    dtOut         = 1.0;  %time interval of output files
    timeSpinUpYr  = 50;   %def: 50. spin-up year  

    % number of amplification of noise (see section 2.2 of Fujita et al. 2024, JGR)
    iterNumMax_org = 10  %default 10
    iterNumMax4Spinup = 50 %default 50

    rngSetting    = 'shuffle' %choose 'shuffle' if you do main simulation, choose 'default' if you want to get the same results for debugging. 
  end

  properties (Access = public)
    status         = 0; 
    currentRun     = 0; 
    currentNode    = 1; 
    currentIterNum = 0; 

    numData        = 1; 
    numCase
    numCase_org 
    numRun         = 1;      
    numNode

    yrInitial      = 1750;  
    yrFinal        = 2015;
    dtOde          = 0.25; %time interval of ode15s 

    iterNumMax

    % model parameter that controls maximum ranges of interannual variations (IAV) of CH4 source and sink parameters
    % def: U(0, 10)% (see section 2.2 of Fujita et al. 2024, JGR for more details)
    % fix: 0.3% (e.g. Egeo, d13Cgeo, dDgeo, wheose time variations are expected to be small)
    paraIAVpercent = struct("def",10,"fix",0.3) 

    % model parameters that control simulations of CH4 source and sink parameters
    controlPara = struct("newDist",'Gauss', ... % type of additional noise distributions: 'Gauss', 'Uniform', 'Cauchy' 
                  "exBoundary",'paraOld') % how to replace the values that exceed the prescribed initial boundary: "paraOld", "NaN"

    useParaTimeVar = true; %true or false
    useDistPrior           %true or false. Set by Test.m/readTargetFile
    useHyperParameter      %true or false. Set by Test.m/readTargetFile

    Info                   %ReadModelParameterInfo
    Spup
    lagFilterTerm

    emsNatr %def "S20": "L07" (Lassey et al., 2007) or "S20" (Saunois et al., 2020) or "dummy" (input zero emission) 
    emsBB   %def "M17": "M17" (van Marle et al. 2017) or "S16" (Schwietzke et al., 2016, 43 Tg/yr constant)
    lossOH 

    issp     = 0;    %def:0 Historical simulation 
    dateLoop = strings(1,2);
    workdir  = "."; 
    randType = "rand";  %"rand", "latin", "Posterior"
    
    inputFile  
    emsCH4File

    atmosOffset = struct('CH4',0,'d13C',0,'dD',0,'D14C',0);

    atmos
    para
    ems
    loss
    judge
    rngList    
    paraTargetList 
    obs 
  end

  properties (SetAccess = private)  
    nameList 
    atmosTargetList  
    nameParaList 
    numParaList  
  end

  methods
% e.g., 1. To run model with prescribed target file: Model(targetDir).main 
%       2. To load existed posterior simulation results: model = Model(targetDir,"Posterior",9)
%       3. To run model using using existed posterior parameters: Model(targetDir,"Posterior",9).main
    function  obj = Model(varargin)  
      addpath constants emission parameter	

      if nargin > 0 
        for i = 1:min(nargin,3)
          if isempty(varargin{i}), continue; end 

          switch i
            case 1, obj.workdir  = varargin{i};
            case 2, obj.randType = varargin{i};
            case 3, obj.status   = varargin{i};
          end
        end
      end

      obj = obj.setup(obj.workdir); 

      if nargin > 3 && mod(nargin,2) == 1
        for k = 4:2:nargin
          obj.(varargin{k}) = varargin{k+1};
        end
      elseif nargin > 3
        error("Please check the nargin is correct!")
      end

      obj.dateLoop(1) = datetimeNow;
    end

%% Setup key model parameter through executing ReadModelParameter and ReadTargetFile.
    function obj = setup(obj,workdir)
      obj.workdir = workdir; 
      obj.Info    = ReadModelParameter(obj.workdir);
      obj.numCase = obj.numCase_org;
      obj.iterNumMax = obj.iterNumMax_org;

      [obj.nameList, obj.atmosTargetList, obj.paraTargetList,...
        obj.numNode, obj.useDistPrior, obj.useHyperParameter,...
        obj.nameParaList,obj.numParaList] ...
        = ReadTargetFile(obj.workdir,obj.inputFile,obj.yrInitial,obj.timeSpinUpYr,obj.atmosOffset).main; 

      assert(obj.atmosTargetList(end).tspan(end) <= obj.yrFinal) 

      obj.writeAtmosTargetListToText; 
      obj.Spup = Spinup(obj.timeSpinUpYr,obj.randType); 
    end

%% Main function of 1-box model.
    function obj = main(obj,varargin) 
      if ~isempty(varargin), obj.currentRun = varargin{1}; end 

      if ismember(obj.randType,["Prior","Posterior","PosteriorMean"]) 
        obj = obj.mainWithoutFiltering; 
        return; 
      end 

      obj = obj.setupMain;       
      for iNode = 2:obj.numNode 
        obj = obj.DoLoopINode(iNode); 
        if obj.status ~= 0, break ; end 
      end 

      %plotTmpFigInModel(obj); %make figure for checking while debugging
      obj = obj.ClearMemorySelect;
      
      if obj.currentRun ~= 0
        save(obj.tmpMatFileName(obj.currentRun),'obj'); 
      end

      disp('End simulation')
    end

%% Setup of main function.
    function obj = setupMain(obj)
      if ~obj.isPosterior, close all; end 
      if obj.currentRun == 0, delete(strcat(obj.workdir,"/tmp*.mat")); end 
      if exist(Smoothing.filenmPostFull(obj,"atmos"),'file') && obj.status == 0, delete(Smoothing.filenmPostFull(obj,"atmos")); end 
      if exist(Smoothing.filenmPostFull(obj,"para"),'file') && obj.status == 0, delete(Smoothing.filenmPostFull(obj,"para")); end   
      obj = obj.InitializeObjectAll; 
    end

%% Initialize main class objects in setupMain & updateIteration 
    function obj = InitializeObjectAll(obj)
      obj.numData  = 1; 
      obj.numCase  = obj.numCase_org; 
      obj.rngList  = RngGenerator(obj); 
      obj.atmos    = Atmosphere(obj); 
      obj.para     = ParameterSourceSink(obj);
      obj.ems      = Emission(obj);
      obj.loss     = Loss(obj);
      obj.judge    = JudgeAtmosToAtmosTarget(obj);

      %if obj.currentNode == 1, writeFilteredResults(obj,obj.currentNode); end 
    end  

%% Sub-function in Loop of main: 
    function obj = DoLoopINode(obj,iNode)
      if isAlreadyReachedYrFinal(obj,iNode), return; end 

      obj.status = 2; 
      if obj.currentIterNum == 0, obj.dateLoop(1) = datetimeNow; end  
      while all(obj.status ~= [0,1]) 
        % Set Initial values
        obj.currentIterNum = obj.currentIterNum + 1; 
        obj                = obj.SetCurrentNodeDataCase(iNode,obj.para(iNode).tspan); 
        obj.rngList(iNode) = obj.SetRng;
        obj.atmos(iNode)   = obj.atmos(iNode).SetAtmosphere(obj,obj.atmos(iNode-1)); 
        obj.para(iNode)    = obj.para(iNode).SetParameter(obj,obj.para(iNode-1)); 
        obj.ems(iNode)     = obj.ems(iNode).SetEmission(obj,obj.para(iNode));
        obj.loss(iNode)    = obj.loss(iNode).SetLoss(obj,obj.para(iNode));
        if iNode == 2 && obj.timeSpinUpYr == 0, obj.status = 0; break; end 

        % Compute atmos
        for gasID = obj.nameList.atmosConc
          [obj.atmos(iNode).tspan,obj.atmos(iNode).(gasID)] = obj.Ode15s(gasID); 
        end
        obj.atmos(iNode) = obj.atmos(iNode).ConvertToIsotopeRatio;

        % Resampling
        if obj.currentIterNum == obj.iterNumMax 
          obj = obj.updateIterationStatus; 
          obj.currentIterNum = 0;
        end
        obj = obj.checkStatus;
      end

      %writeFilteredResults(obj,iNode); 
      obj.dateLoop(2) = datetimeNow;

      function y = isAlreadyReachedYrFinal(obj,iNode)
        y = iNode == obj.numNode && obj.atmos(iNode-1).tspan(end) == obj.yrFinal;
      end
    end   

%% Setup of numData & numCase   
    function obj = SetCurrentNodeDataCase(obj,iNode,tspan) 
      obj.numData = length(tspan);    
      [obj.numCase,obj.currentNode] = obj.updateNumCaseAndCurrentNode(iNode); 
      
      if obj.currentIterNum == 1, disp(strcat("iNode = ",num2str(iNode))); end
      if iNode == obj.numNode
        disp(strcat("Integral period: -",num2str(obj.yrFinal)));        
      else
        disp(strcat("Integral period: ",strrep(num2str(obj.atmosTargetList(iNode-1).tspan),"  ","-")));
      end

      if obj.currentNode == 2 && obj.Spup.step == 0 && ~isempty(obj.Spup.type) 
        %obj.Info.iterNumMax = obj.iterNumMax4Spinup; 
        obj.iterNumMax = obj.iterNumMax4Spinup; 
      end
    end

%% Sub-function of SetCurrentNodeDataCase
    function [numCase,currentNode] = updateNumCaseAndCurrentNode(obj,iNode) 
      if obj.numParaList(iNode).fDynamic > 0  
        numCase = obj.numCase_org;
      else
        numCase = obj.judge(iNode-1).numMatch.all;
      end
      currentNode = iNode;
    end   

%% Update obj.status for controlling loop
    function obj = updateIterationStatus(obj) 
      iNode = obj.currentNode;
      disp(strcat("iWhile = ",num2str(obj.currentIterNum)," for iteration")); 

      obj = obj.updateResultsResIdx; 

      if obj.Spup.needsCalcAtmosEquilibrium && obj.judge(iNode).numMatch.all < 3
        if ~obj.isOverTimeFromLastStep(5)
          obj = obj.InitializeObjectAll; 
          delete(obj.tmpIterObjMatFileNameInLoop("*"));
          disp('Particle number during first Spinup < 3, then retry prediction!');
          return
        else
          error(strcat("Cannot find the solution at iNode = ",num2str(iNode))); 
        end
      end

      if iNode == 2 && obj.Spup.step == 0 && ~isempty(obj.Spup.type)
        obj.Spup = obj.Spup.saveMatFileTmpParaTarget(obj); 

        if ismember(obj.Spup.type,["noise","Cauchy"])
          %obj.Info.iterNumMax = min([max([obj.numCase,obj.iterNumMax4Spinup]),50]); 
          obj.iterNumMax = min([max([obj.numCase,obj.iterNumMax4Spinup]),50]); 
        end

        obj.numCase = obj.numCase_org; 
        obj.Spup.step = 1;
        disp('Complete first spin-up.');

      elseif iNode == 2 && obj.Spup.step == 1
        if ismember(obj.Spup.type,["noise","Cauchy"])
          obj.iterNumMax = obj.iterNumMax_org; 
          obj.Spup.needsCalcAtmosEquilibrium = false;
        end

        obj.status = 0;
        obj.Spup.step = 0;
        disp('Complete second spin-up. Proper initial values have been optimized.');

      else
        obj.status = 0;
        disp('Update object from amplified predicted distribution. Iteration has been completed.');
      end
    end

    function y = isOverTimeFromLastStep(obj,valMinutes) 
      y = duration(datetime(datetimeNow)-datetime(obj.dateLoop(2))) > duration(0,valMinutes,0);
    end

%% Update each object from observation likelihood at iNode
    function obj = updateResultsResIdx(obj)
      iNode = obj.currentNode;
      if iNode == obj.numNode, return; end 

      % obj.judge(iNode).resIdx is created based on the likelihood.
      % The initial ensembles of atmos(iNode+1) & para(iNode+1) is created as atmos(iNode)(end,obj.judge(iNode).resIdx) in the next loop.
      [atmosComb,paraComb]  = obj.combineIterAtmosPara;
      disp('Load tmp_atmos & tmp_para')
      obj.judge(iNode) = obj.judge(iNode).getJudgeMatchAtmosTarget(obj,atmosComb,paraComb); 
      disp(obj.judge(iNode).numMatch);
      obj.judge(iNode).parentIdx = obj.judge(iNode).convertJudgeIsMatchAllToParentIdx(obj); 
      obj.judge(iNode) = obj.judge(iNode).updateResIdx; 

      obj.numCase      = obj.judge(iNode).numMatch.all; 
      assert(obj.numCase > 0, strcat("Error: iNode = ",num2str(iNode), ", numCase = 0"))
      obj.atmos(iNode) = obj.judge(iNode).updateResults(obj.judge(iNode).isMatch.all,atmosComb);
      obj.para(iNode)  = obj.judge(iNode).updateResults(obj.judge(iNode).isMatch.all,paraComb);
      obj.ems(iNode)   = obj.ems(iNode).SetEmission(obj,obj.para(iNode));
      obj.loss(iNode)  = obj.loss(iNode).SetLoss(obj,obj.para(iNode));

      % if obj.currentRun == 0 && obj.Spup.step == 1 
      %   plotFig(obj,paraComb);
      % end

      function plotFig(obj,paraComb)
        fig = figure; t = tiledlayout(4,5); 
        for name = obj.nameList.para, nexttile;
          [paraMin,paraMax]=MyFunc.SetParaMinMax(name);
          edgePara = paraMin:(paraMax-paraMin)/20:paraMax;
          histogram(paraComb.(name)(end,:),edgePara,'Normalization','probability'); hold on  %proposed distritbuion
          histogram(obj.para(2).(name)(end,obj.judge(2).resIdx),edgePara,'Normalization','probability');  %filtered distritbuion
          ylabel(MyFunc.nameLabel(name,"para"));
        end
        t.Title.String = {strcat("Spinup noise prercent = ",num2str(obj.Spup.paraNoisePercent*100),"%"); ...
                         strcat("iterNumMax: ",num2str(obj.iterNumMax),", numCase: ",num2str(obj.judge(2).numMatch.all),"/",num2str(obj.numCase_org))};
        
        fig.Position(3) = 900;
        saveas(fig,strcat(obj.workdir,"/HistSpup_",obj.Spup.type,num2str(obj.Spup.paraNoisePercent*100),"Per.png"))
        disp("")
      end

    end

%% Combine temporally matfile created during loop in main  
    function [atmos,para] = combineIterAtmosPara(obj)
      atmos = obj.atmos(obj.currentNode);
      para  = obj.para(obj.currentNode);
      for iterNum = obj.iterNumMax-1:-1:1
        filename = obj.tmpIterObjMatFileNameInLoop(iterNum); 
        tmp = load(filename);
        atmos = combine(tmp.atmos, atmos, unique([obj.nameList.atmos, obj.nameList.atmosConc]));
        para  = combine(tmp.para, para, obj.nameList.para); 
      end

      function obj = combine(tmpObj, obj, nameList, varargin)
        for name = nameList
          obj.(name) = horzcat(tmpObj.(name), obj.(name)); 
          if isa(obj,'ParameterSourceSink') && obj.ismemberHyper(name) 
            obj.hyper.(name) = horzcat(tmpObj.hyper.(name), obj.hyper.(name));
          end     
          if isa(obj,'ParameterSourceSink') && obj.containsIAVpercentEach
            obj.IAVpercent.(name) = horzcat(tmpObj.IAVpercent.(name), obj.IAVpercent.(name));
          end        
        end
      end

    end

%% 
    function obj = checkStatus(obj)
      switch obj.status

        % Normal termination: more than one particles remained and the number of spinups (obj.iterNumMax) is reached 
        case 0 
          disp(strcat(datetimeNow,' iteration has been completed.')); 
        
        % Error termination: no particle remained 
        case 1 
          disp(strcat(datetimeNow,' there is no particle that match the observational target under ', num2str(obj.iterNumMax) ,'  iteration. Please change the model setting.'));        

        % Still in the process of finding particles that match observations
        case 2 
          obj.saveAtmosParaTmpMat; 
          if obj.currentNode == 2 
            obj = obj.InitializeObjectAll;
          end
      end
    end

%% Save tmporally mat file during iteration resampling
    function saveAtmosParaTmpMat(obj)
      atmos = obj.atmos(obj.currentNode); 
      para = obj.para(obj.currentNode);
      save(obj.tmpIterObjMatFileNameInLoop(obj.currentIterNum),'atmos','para'); clear atmos para
    end

%%
    function y = tmpIterObjMatFileNameInLoop(obj,iterNum)
      fileHead = strcat(obj.workdir,'/tmp_IterObj_',num2str(obj.currentRun,'%05.0f'),'_Loop');
      if isa(iterNum,'double')
        y = strcat(fileHead,num2str(iterNum,'%02.0f'),'.mat');
      elseif isa(iterNum,'string')
        y = strcat(fileHead,iterNum,'.mat');
      else
        error("Wrong input format: iterNum")
      end
    end

%% Perform Ode15s without resampling. 
%  obj.randType = 'Prior' or 'Posterior' or 'PosteriorMean'
%  e.g. model = Model('.','PosteriorMean').main;
    function obj = mainWithoutFiltering(obj)
      obj = obj.SetupWithoutFiltering; 

      dtOde_org = obj.dtOde;
      obj.dtOde = 1/12; %required to properly integrate with ode15s over the period when the emission changes are rapid year by year
      for gasID = obj.nameList.atmosConc
        [obj.atmos.tspan,obj.atmos.(gasID)] = obj.Ode15s(gasID); 
      end
      obj.dtOde = dtOde_org; 

      obj.atmos = obj.atmos.ConvertToIsotopeRatio;
      obj = obj.ClearMemorySelect;
      save(obj.tmpMatFileName(obj.currentRun),'obj');
      %plotFig(obj);

      function plotFig(obj)
        obj.atmos.CH4 = obj.atmos.CH4./MassConstants.factor;
        figure; tiledlayout(2,2); for gas = obj.nameList.atmos, nexttile; plot(obj.atmos.tspan,obj.atmos.(gas)); end
        figure; tiledlayout(4,5); 
        for paraID = obj.nameList.para, nexttile;
          [paraMin,paraMax]=MyFunc.SetParaMinMax(paraID);
          histogram(obj.para.(paraID),paraMin:(paraMax-paraMin)/10:paraMax);
        end
      end
    end

%% Setup of mainWithoutFiltering
    function obj = SetupWithoutFiltering(obj,varargin) 
      obj.status = 10; 
      switch obj.randType
        case {"Prior"} 
          obj.numNode = 1; 
          obj.atmos   = Atmosphere(obj);  
        case {"PosteriorMean"} 
          obj.numNode = 1; obj.numCase = 1; 
          obj.atmos   = Atmosphere(obj);  
        case {"Posterior"}
          obj.atmos = Atmosphere(obj); 
          obj.atmos = obj.atmos.SetAtmosphereInitial(obj);
      end

      obj.numData = length(obj.atmos.tspan)+1;
      obj.para    = ParameterSourceSink(obj); 
      obj.ems     = Emission(obj,varargin{:});
      obj.loss    = Loss(obj);
    end 

%% Create random seed object of data storage for each iNode
    function rngList = RngGenerator(obj)
      rng(obj.rngSetting);               %Set current seed
      rngList              = rng;        %Allocate current seed
      rngList(obj.numNode) = rngList(1); %Pre-allocation of class arrays
    end

%% Set or call random seed in obj.rngList
    function y = SetRng(obj)
      if isempty(obj.rngList(obj.currentNode).Type) && obj.currentNode <= 2 
        y = obj.rngList(1);   
      elseif isempty(obj.rngList(obj.currentNode).Type) && obj.currentNode >  2
        rng(obj.rngSetting); 
        y = rng;             
      else
        y = obj.rngList(obj.currentNode);
      end
    end

%% Clear variables that will not be used in the analysis or that can be reconstructed afterward
    function obj = ClearMemorySelect(obj)
      obj.ems   = []; 
      obj.loss  = []; 
    end

%% Load smoothed file objects_SmoothSummary_Full.mat or objects_SmoothSummaryPostFull_*.mat (if exist)
    function obj = loadPosteriorEnsemble(obj,varargin)
      obj.randType = "Posterior"; 
      obj.status   = 9;
      obj.atmos    = Atmosphere(obj,'ConvertToIsotopeMoleFraction','DevideMassFactorToMoleFraction'); 
      obj.para     = ParameterSourceSink(obj);
      obj.ems      = Emission(obj,varargin{:});
      obj.loss     = Loss(obj);   
    end

%% Load objects_*Summary.mat (filtered distributions for each iRun and iNode)
    function obj = loadFilteredData(obj,propNameList,varargin)
      disp(strcat("Start loadFilteredData ",datetimeNow))

      for propName = [propNameList,"judge"] %judge is always required to extract smoothed data from objects_*Summary.mat
        % To avoid twice reading 
        if strcmp(propName,"paraIAV")
          propName = "para";
        elseif ~isempty(obj.(propName))
          disp(strcat('model.',propName,' is already allocated. Thus continue.')); 
          continue; 
        end 

        filename = strcat(obj.workdir,'/objects_',propName,'Summary.mat'); 

        switch propName
          case {'ems','loss'}
            obj = obj.reconstructObject(propName);
          otherwise
            obj = obj.loadMatFile(propName,filename); 
        end
  
        % Load selected iRunIndex target (e.g., 1:5 of 1:10)       
        if ~isempty(varargin)
          iRunIdxSelect  = varargin{1}; 
          obj.(propName) = obj.(propName)(iRunIdxSelect,:); 
        end 
      end
      disp(strcat("End loadFilteredData ",datetimeNow))
    end

%% Create ems or loss object from objects_paraSummary.mat
    function obj = reconstructObject(obj,propName)
      disp(strcat("Start reconstructObject ",propName," ",datetimeNow))
      obj.status = 9; 
      obj = obj.loadFilteredData("para");
      switch propName
        case 'ems'
          obj.(propName) = Emission(obj); 
        case 'loss'
          obj.(propName) = Loss(obj); 
      end 
      disp(strcat("Finish reconstructObject ",propName," ",datetimeNow))
    end

%% 
    function obj = loadMatFile(obj,propName,filename)  
      if strcmp(propName,"judge") && ~exist(filename,'file'), return; end 
      disp(strcat("Start loadMatFile ",filename," ",datetimeNow))

      %Replace filename if it is ssp simulation 
      if obj.issp ~=0 && strlength(char(obj.randType)) == 4 && strcmp(obj.randType,"OSSE")
        filename = strrep(filename,".mat",strcat(MyNameList.ssp(obj.issp),"OSSE.mat"));
      elseif obj.issp ~=0
        filename = strrep(filename,".mat",strcat(MyNameList.ssp(obj.issp),".mat"));
      end

      obj.(propName) = matfile(filename).var; 
      disp(strcat("Finish loadMatFile ",filename," ",datetimeNow))
    end

%% 
    function y = tmpMatFileName(obj,iRun,varargin)
      switch obj.randType
        case {"Prior","Posterior","PosteriorMean"}
          header = strcat("tmp_",obj.randType);
        otherwise
          if isempty(varargin), header = 'tmp_objects'; 
          else, header = varargin{1}; 
          end
      end

      y  = strcat(obj.workdir,"/",header,"_",num2str(iRun,'%05.0f'),'.mat');
    end

%%
    function y = isPosterior(obj)
      y = strlength(char(obj.randType)) == 9 && strcmp(obj.randType,"Posterior"); 
    end

%%
    function y = isPosteriorOSSE(obj)
      y = strlength(char(obj.randType)) == 4 && strcmp(obj.randType,"OSSE") && obj.status == 9;
    end

%% 
    function y = isLagSmoother(obj)
      y = obj.status == 8; 
    end

%% 
    function y = isReconstruction(obj) 
      y = obj.status == 9; 
    end

%% 
    function y = isRunWithoutFiltering(obj)
      y = obj.status == 10;
    end

%% 
    function y = isPriorOnly(obj)
      y = obj.isRunWithoutFiltering && strcmp(obj.randType,'Prior');
    end

%% 
    function y = isPosteriorMeanOnly(obj)
      y = obj.isRunWithoutFiltering && strcmp(obj.randType,'PosteriorMean');
    end

%%
    function writeAtmosTargetListToText(obj)
      filename = strcat(obj.workdir,"/atmosTargetList.txt");
      if exist(filename,"file"), return; end
    
      fileID = fopen(filename,'w');
      for i = 1:length(obj.atmosTargetList)
        atmosTar = obj.atmosTargetList(i);
        tmp = [atmosTar.termName,atmosTar.tspan,struct2cell(atmosTar.CH4).',...
          struct2cell(atmosTar.d13C).',struct2cell(atmosTar.dD).',struct2cell(atmosTar.D14C).'];
        fprintf(fileID,'%s\t',fillmissing(tmp,'constant',"NaN")); fprintf(fileID,'\n'); 
      end
      fclose(fileID);
    end

%% 
    function  [namePara,numPara] = extractNameNumParaAllTerm(obj)    
      namePara = struct('termName',"",'all',"",'f0',"",'f1',"",'f2',"",'f3',"",'f4',"",'f5',"",'f6',"",'f7',"",'f8',"",'f9',"",'fDynamic',"",'hyper',"");
      numPara  = struct('termName',"",'all',0,'f0',0,'f1',0,'f2',0,'f3',0,'f4',0,'f5',0,'f6',0,'f7',0,'f8',0,'f9',0,'fDynamic',0,'hyper',0);
      for iNode = 2:obj.numNode, namePara(iNode) = namePara(1); numPara(iNode) = numPara(1); end      

      for iNode = 1:obj.numNode-1 
        namePara(iNode).all  = obj.nameList.para;
        numPara(iNode).all   = length(namePara(iNode).all);

        if iNode == 1         
          paraTarget               = obj.paraTargetList(1); 
          namePara(iNode).termName = 'initial';
          numPara(iNode).termName  = 'initial';         
        else
          paraTarget               = obj.paraTargetList(iNode-1);
          namePara(iNode).termName = paraTarget.termName;
          numPara(iNode).termName  = paraTarget.termName;
        end
        
        for name = namePara(iNode).all
          fX = strcat('f',num2str(paraTarget.(name).flag));
          [namePara(iNode),numPara(iNode)] = setForNumName(fX,name,namePara(iNode),numPara(iNode));
          if paraTarget.(name).flag ~= 0, [namePara(iNode),numPara(iNode)] = setForNumName('fDynamic',name,namePara(iNode),numPara(iNode)); end
          if obj.useHyperParameter && paraTarget.(name).setHyper,  [namePara(iNode),numPara(iNode)] = setForNumName('hyper',name,namePara(iNode),numPara(iNode)); end %2022.5.30
        end   
      end
      namePara(end) = namePara(end-1); namePara(end).termName = 'final';
      numPara(end)  = numPara(end-1);  numPara(end).termName  = 'final';

      function [namePara,numPara] = setForNumName(fX,name,namePara,numPara)
        numPara.(fX) = numPara.(fX) + 1;
        namePara.(fX)(numPara.(fX)) = name;
      end
    end

%% 
    function [t,y] = Ode15s(obj,gasID,varargin) 
      opts  = odeset('MaxStep',obj.dtOde); 
      
      if isempty(varargin) 
        iNode      = obj.currentNode;
        tspanAtmos = obj.atmos(iNode).tspan;
        tspanEms   = obj.ems(iNode).tspan;
        emsOde     = obj.ems(iNode);
        lossOde    = obj.loss(iNode);
        y0         = obj.atmos(iNode).(gasID)(1,:).';    
        cname      = "tot"; 
      else
        tspanAtmos = varargin{1}; 
        atmosOde   = varargin{2}; 
        tspanEms   = varargin{3}; 
        emsOde     = varargin{4};    
        lossOde    = varargin{5}; 
        cname      = varargin{6}; 
        y0         = atmosOde.(gasID)(1,:).'; 
        clear varargin
      end

      validIdx = ~any(isnan(emsOde.(gasID).(cname))) & ~any(isnan(lossOde.(gasID).tot));
      if sum(validIdx) ~= obj.numCase
        emsOde.(gasID).(cname)  = emsOde.(gasID).(cname)(:,validIdx);
        lossOde.(gasID).(cname) = lossOde.(gasID).(cname)(:,validIdx);
        y0 = y0(validIdx); %If at least one NaN is included in the emsOde and/or lossOde rows, do not calculate the integral for that index.
      end

      if obj.Spup.needsCalcAtmosEquilibrium 
        [t,y] = obj.atmosEquilibrium(gasID,emsOde,lossOde,tspanAtmos,y0); 
      else
        [t,y] = ode15s(@(t,y) obj.massBalance(gasID,tspanEms,emsOde,lossOde,t,y,cname), tspanAtmos, y0, opts); t = t.'; 
      end

      if sum(validIdx) ~= obj.numCase 
        y_tmp = y; 
        y = NaN(length(t),obj.numCase);
        y(:,validIdx) = y_tmp; %keep the total array size as obj.numCase
      end
    end

%% 
    function y = get.numCase_org(obj)
      y = obj.Info.numCase;
    end

%%
    function y = get.emsCH4File(obj) 
      y = obj.Info.inventory + ".txt";
    end
%%
    function y = get.inputFile(obj)
      y = obj.Info.atmosTargetFile;
    end

%%
    function y = get.emsNatr(obj) 
      y = obj.Info.emsNatr;
    end

%%
    function y = get.emsBB(obj) 
      y = obj.Info.emsBB;
    end

%% 
    function y = get.lossOH(obj) 
      y = obj.Info.lossOH;
    end

  end

%%
  methods(Static)
    function [t,y] = atmosEquilibrium(gasID,ems,loss,tspanAtmos,y0) 
      y      = zeros(2,length(y0));
      y(1,:) = y0;
  
      iLast  = size(ems.(gasID).tot,1);
      y(2,:) = ems.(gasID).tot(iLast,:)./loss.(gasID).tot(iLast,:);
      t      = [tspanAtmos(1),tspanAtmos(end)];  
    end

    %%
    function y = massBalance(gasID,tspanEms,ems,loss,t,y,cname)
      y = interp1(tspanEms, ems.(gasID).(cname) - y.'.*loss.(gasID).tot,  t, 'linear').';
    end

  end

end

%%
function y = datetimeNow()
  y = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
end
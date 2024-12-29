% Contains particle filter code for filtering. It compares
% ensembles of simulated atmospheric CH4, d13C, dD, and 
% D14C with their observational targets and then performs 
% resampling using the likelihood function of the observations.
% Ryo Fujita 2024

classdef JudgeAtmosToAtmosTarget
  properties
    termName
    isMatch
    numMatch
    likelihood  
    resIdx      
    parentIdx   
  end

%% 
  methods
    function obj = JudgeAtmosToAtmosTarget(model)
    if nargin ~= 0
      if model.isRunWithoutFiltering 
        numSize = 1; 
      else
        numSize = size(model.atmosTargetList,2)+2;
      end

      numCase = model.numCase;

      for i = 1:numSize
        switch i
          case 1
            obj(i).termName = "initial";
          case numSize
            obj(i).termName = "final";
          otherwise
            obj(i).termName = model.nameList.term(i-1);
        end

        obj(i).isMatch = struct("all",[],"CH4",[],"d13C",[],"dD",[],"D14C",[]);
        obj(i).numMatch = struct("all",0,"CH4",0,"d13C",0,"dD",0,"D14C",0);
        obj(i).likelihood = struct("all",[],"CH4",[],"d13C",[],"dD",[],"D14C",[]);
        obj(i).resIdx = 1:numCase;
        obj(i).parentIdx = 1:numCase;
      end
      
      for gasID = ["CH4","d13C","dD","D14C"] 
        obj(1).isMatch.(gasID)   = true(1,numCase); 
        obj(1).numMatch.(gasID)  = numCase;
      end
      obj(1).isMatch.all   = true(1,numCase);
      obj(1).numMatch.all  = sum(obj(1).isMatch.all,2); 
    end
    end

%% 
    function obj = getJudgeMatchAtmosTarget(obj,model,varargin)
      if model.currentNode == model.numNode 
        obj = obj.getJudgeAtmosphere_iNodeLast(model.atmos(model.currentNode));
        return
      end
 
      if ~isempty(varargin) && length(varargin) == 2
        atmos = varargin{1};
        para  = varargin{2};
        clear varargin
      else
        atmos = model.atmos(model.currentNode);
        para  = model.para(model.currentNode);
      end

      for gasID = ["CH4","d13C","dD","D14C"] 
        if atmos.target.(gasID).flag == 0
          atmos = obj.convertAtmosTargetFlag0(atmos,gasID); 
        end

        obj.isMatch.(gasID)  = obj.getIsMatch(atmos,gasID);
        obj.numMatch.(gasID) = sum(obj.isMatch.(gasID));
        obj = obj.computeLikelihood(model,atmos,gasID,para);
      end

      obj.likelihood.all = obj.likelihood.CH4 + obj.likelihood.d13C + ...
        obj.likelihood.dD + obj.likelihood.D14C;  
      obj.resIdx         = obj.resampleParticle('all',model.numCase_org);
      obj.isMatch.all    = obj.convertResIdxToIsMatchAll; 
      obj.numMatch.all   = sum(obj.isMatch.all);
    end

%% 
    function obj = getJudgeAtmosphere_iNodeLast(obj,atmos)
      for gasID = ["CH4","d13C","dD","D14C"]
        obj.isMatch(1).(gasID)   = true(1,size(atmos.(gasID),2)); 
        obj.numMatch(1).(gasID)  = sum(obj.isMatch(1).(gasID),2); 
      end
      obj.isMatch(1).all  = obj.isMatch(1).CH4 == 1 & obj.isMatch(1).d13C == 1 & ...
        obj.isMatch(1).dD == 1 & obj.isMatch(1).D14C == 1;
      obj.numMatch(1).all = sum(obj.isMatch(1).all);
    end

%% 
    function obj = computeLikelihood(obj,model,atmos,gasID,para)
      if strcmp(gasID,'CH4'), atmos.(gasID) = atmos.(gasID)./MassConstants.factor; end
      obj.likelihood.(gasID) = zeros(1,size(atmos.(gasID),2));

      if atmos.target.(gasID).flag == 2 ...
        && sum(ismember(fieldnames(atmos.target.(gasID)),["ave","sdev"])) == 2 
        A = -log(sqrt(2 * pi) * atmos.target.(gasID).sdev); 
        B = - 0.5 / (atmos.target.(gasID).sdev.^2);
        
        logic = para.isWithinPriorRangeIAV(model);
        for k = 1:size(atmos.(gasID),2)
          if logic(k) 
            % Square of the difference between predicted and observed values
            D2 = (atmos.(gasID)(end,k) - atmos.target.(gasID).ave)' * ...
              (atmos.(gasID)(end,k) - atmos.target.(gasID).ave); 
            
            % log likelihood of particle (k)
            obj.likelihood.(gasID)(k) =  A + B * D2;
          else 
            % Likelihood is minus infinity (a value that does not match the target at all)
            obj.likelihood.(gasID)(k) = -Inf;
          end
        end

      else
        logic = obj.isMatch.(gasID) & para.isWithinPriorRangeIAV(model);
        for k = 1:size(atmos.(gasID),2) 
          if logic(k)
            obj.likelihood.(gasID)(k) = 0;
          else 
            obj.likelihood.(gasID)(k) = -Inf;
          end
        end
      end
    end

%% 
    function resIdx = resampleParticle(obj,gasID,numCase_org)
      % Calculating Cumulative Distribution
      L = exp(obj.likelihood.(gasID) - max(obj.likelihood.(gasID)));
      Q = L / sum(L, 2); 
      R = cumsum(Q, 2);  

      if ~any(isfinite(R))
        disp('All likelihood are infinit! Need resampling.');
        resIdx = []; 
        return
      end  

      numCase = numCase_org; 

      % Assign the number of samples according to the likelihood ratio (integer portion of Q*numCase)
      m = floor(Q.*numCase); 
      idxList    = find(m~=0);
      numIdxList = m(idxList);

      resIdx = zeros(1,numCase); ist=1;
      for iList = 1:length(idxList)
        iend = ist + numIdxList(iList) -1;
        resIdx(ist:iend) = idxList(iList); 
        ist = iend + 1; 
      end
      
      % Random sampling for the remaining number of samples
      if sum(m) ~= numCase
        % The minority portion of the original weights was used as the new weights
        Q_residual = (Q.*numCase - m)./(numCase-sum(m)); 
        R_redisual = cumsum(Q_residual, 2);
        [~, I] = histc(rand(1, numCase-sum(m)), R_redisual);
        resIdx(ist:numCase) = I + 1;
      end
    end

%% Output the index number specified in created resIdx as true
    function y = convertResIdxToIsMatchAll(obj)
      y = false(1,size(obj.likelihood.all,2));
      if isempty(obj.resIdx)
        return; 
      end

      for idx = unique(obj.resIdx)
        y(idx) = true;
      end
    end

%% 
    function parentIdx = convertJudgeIsMatchAllToParentIdx(obj,model)
      parentIdx = mod(find(obj.isMatch.all == 1),model.numCase_org); 
      parentIdx(parentIdx==0) = model.numCase_org;
    end

%% 
    function y = isMatchExclCH4(obj)
      y = obj.isMatch.d13C == 1 & obj.isMatch.dD == 1 & obj.isMatch.D14C == 1;
    end

%% 
    function obj = updateResIdx(obj)
      idxNew = 0;
      for idx = unique(obj.resIdx)
        idxNew = idxNew + 1;
        obj.resIdx(obj.resIdx == idx) = idxNew;
      end
    end

  end

%% 
  methods (Static)
%%
    function inout = updateResults(judgeLogic,inout,varargin) 
      switch class(inout)
        case 'Atmosphere'
          fieldList = ["CH4","d13C","dD","D14C","C13","CHD","C14"];
        case {'ParameterSourceSink','ParaTimeVar','HyperParameter'} 
          if ~isempty(varargin)
            fieldList = varargin{1};
          else
            listTmp   = string(fieldnames(ParameterSourceSink)).';
            fieldList = listTmp(~ismember(listTmp,["tspan" "termName" "target" "IAVpercent" "hyper"])); 
          end
      end

      switch class(inout)
        case {'Atmosphere','ParameterSourceSink','ParaTimeVar','HyperParameter'}
          for name = fieldList
            if isempty(inout.(name))
              continue; 
            end
 
            if isempty(judgeLogic)
              inout.(name) = []; 
              continue; 
            end
 
            inout.(name) = inout.(name)(:,judgeLogic); 

            if isa(inout,'ParameterSourceSink') && inout.ismemberHyper(name) 
              inout.hyper = JudgeAtmosToAtmosTarget.updateResults(judgeLogic,inout.hyper,name);
            end
          end

          if isa(inout,'ParameterSourceSink') && inout.containsIAVpercentEach 
            inout.IAVpercent = JudgeAtmosToAtmosTarget.updateResults(judgeLogic,inout.IAVpercent);
          end
        case {'Emission'}
          for gasID = ["CH4","C13","CHD","C14"]
            fieldList = string(fieldnames(inout.(gasID))).';
            for name = fieldList
              if isempty(judgeLogic)
                inout.(gasID).(name) = []; 
                continue; 
              end 
              inout.(gasID).(name) = inout.(gasID).(name)(:,judgeLogic); 
            end
          end
        case {'Loss'}
          for gasID = ["CH4","C13","CHD","C14"]
            fieldList = string(fieldnames(inout.(gasID))).';
            for name = fieldList
              if isempty(judgeLogic)
                inout.(gasID).(name) = []; 
                continue;   
              end
              inout.(gasID).(name) = inout.(gasID).(name)(:,judgeLogic); 
            end
          end
      end
    end

%% atmospheric targets that not used for optimization 
% (set as flag = 0 in the target file)
    function atmos = convertAtmosTargetFlag0(atmos,gasID)
      switch gasID
        case 'CH4'
          if atmos.target.(gasID).flag == 0
            atmos.target.(gasID).min = 0; 
            atmos.target.(gasID).max = 9999; 
          end
        otherwise
          if atmos.target.(gasID).flag == 0
            atmos.target.(gasID).min = -1000; 
            atmos.target.(gasID).max = 1000; 
          end
      end
    end

%% 
    function y = getIsMatch(atmos,gasID)
      switch gasID
        case 'CH4'
          y = atmos.target.CH4.min <= atmos.CH4(size(atmos.CH4,1),:)./MassConstants.factor & ...
            atmos.CH4(size(atmos.CH4,1),:)./MassConstants.factor  <= atmos.target.CH4.max; 
        otherwise
          y = atmos.target.(gasID).min <= atmos.(gasID)(size(atmos.(gasID),1),:) & ...
            atmos.(gasID)(size(atmos.(gasID),1),:)  <= atmos.target.(gasID).max;
      end
    end

%%
    function judgeComb = combineIsMatch(nameID,judge)
      numRun  = size(judge,1);
      numCase_org = 0; 

      for iRun = 1:numRun
        numCase_org = numCase_org + size(judge(iRun).isMatch.(nameID),2); 
      end

      judgeComb = zeros(1,numCase_org);

      iBegin = 1; 
      for iRun = 1:numRun
        iEnd = iBegin + size(judge(iRun).isMatch.(nameID),2) - 1;
        judgeComb(:,iBegin:iEnd) = judge(iRun).isMatch.(nameID);
        iBegin = iEnd + 1;
      end
    end

%% 
    function y = shortName()
      y = "judge";
    end

  end % End methods(Static)
end
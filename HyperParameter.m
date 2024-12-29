% Handles the hyperparameters.
% Ryo Fujita 2024

classdef HyperParameter < ParameterSourceSink
  %flag
  % 1: extract from a new prior distribution by each time step
  % 2: update parameters from one period ago
  % 3: update parameters from one period ago with small noise (default)
  properties
    flag  = 3;
    noise = 0.01; 
  end

  methods
    function obj = HyperParameter(model)
      if nargin == 0
        return; 
      end
 
      obj(model.numNode) = obj;
      for iNode = 1:model.numNode
        model.currentNode = iNode;
        switch iNode
          case 1
            obj(iNode).termName = "initial";
            obj(iNode).tspan    = model.paraTargetList(1).tspan(1);
            obj(iNode).target   = model.paraTargetList(1);
          case 2
            obj(iNode).termName = "sp";
            obj(iNode).tspan    = model.paraTargetList(1).tspan(2);
            obj(iNode).target   = model.paraTargetList(1);
            obj(iNode)          = obj(iNode).Initialize(model);
          case model.numNode
            obj(iNode).termName = "final";
          otherwise
            obj(iNode).termName = model.nameList.term(iNode-1);
            obj(iNode).tspan    = model.paraTargetList(iNode-1).tspan(2);
            obj(iNode).target   = model.paraTargetList(iNode-1);
        end
      end
    end

%% 
    function para = setup(obj,para)
      for iNode = 1:size(obj,2)
        para(iNode).hyper = obj(iNode);
      end
    end

%%
    function obj = SetParameter(obj,model)
      if model.currentNode == 1
        % Create initial distrubutions under iNode = 1
        obj = obj.Initialize(model); 
        return
      end

      switch obj.flag
        %case 0, obj = obj.Fix(model); % Use constant value IAVpercent (TBD)
        case 1, obj = obj.Initialize(model); 
        case 2, obj = obj.Update(model);     
        case 3, obj = obj.UpdateInovation(model); 
      end
    end

%%
    function obj = UpdateInovation(obj,model) 
      iNode  = model.currentNode;
      resIdx = model.judge(iNode-1).resIdx; 
      rng(model.rngSetting); 

      for name = model.nameList.para
        switch obj.target.(name).setHyper
          case 0
            obj.(name) = [];
          case 1
            obj_old = model.para(iNode-1).hyper.(name)(1,resIdx);
            obj.(name) = obj_old + obj.noise .* (-1 + 2*rand(1,length(resIdx)));

            logicBoundary = obj.(name) < obj.target.(name).p1 | obj.(name) > obj.target.(name).p2;
            obj.(name)(logicBoundary) = obj_old(logicBoundary); 
        end
      end
    end

%% 
    function obj = Initialize(obj,model) 
      numPara = model.numParaList(model.currentNode); 
      rng(model.rngSetting); 
      randArray = obj.CreateRandArray(model.numCase,numPara.hyper,"latin"); 

      iRand = 1;
      for name = model.nameList.para
        switch obj.target.(name).setHyper
          case 0 
            obj.(name) = []; 
          case 1
            obj.(name) = obj.ComputeParaFullRange(name,randArray(iRand,:));
            iRand = iRand + 1;
          otherwise
            error('Not yet prepared!')
        end 
      end
    end

%%
    function y = ComputeParaFullRange(obj,name,randArray)
      y = obj.target.(name).p1 + (obj.target.(name).p2-obj.target.(name).p1) .* randArray; 
    end
  end
end

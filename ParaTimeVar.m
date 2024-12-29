% Handles the time variations of parameters (i.e. hyperparameter).
% Ryo Fujita 2024

classdef ParaTimeVar < ParameterSourceSink
  %flag
  % 1: extract from a new prior distribution by each time step
  % 2: update parameters from one period ago
  % 3: update parameters from one period ago with small noise (default)
  properties
    flag  = 3; 
    noise = 1; %percent 
  end

  methods
    function obj = ParaTimeVar(model)
      if nargin == 0, return; end 

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
    function para = SetIAVpercent(obj,para)
      for iNode = 1:size(obj,2)
        para(iNode).IAVpercent = obj(iNode);
      end
    end

%% 
    function obj = SetParameter(obj,model)
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

      if iNode == 3
        obj = obj.Initialize(model); 
        return; 
      end

      resIdx = model.judge(iNode-1).resIdx; 
      rng(model.rngSetting);

      for name = model.nameList.para
        switch obj.target.(name).flag
          case {0,1,2}
            obj.(name) = []; 
          case 4
            obj_old = model.para(iNode-1).IAVpercent.(name)(1,resIdx);
            obj.(name) = obj_old - obj.noise .* (-1 + 2*rand(1,length(resIdx)));

            logicBoundary = obj.(name) < 0 | obj.(name) > model.paraIAVpercent.def;
            obj.(name)(logicBoundary) = obj_old(logicBoundary);
        end
      end
    end

%% 
    function obj = Update(obj,model)
      iNode  = model.currentNode;
      resIdx = model.judge(iNode-1).resIdx; 
      for name = model.nameList.para
        switch obj.target.(name).flag
          case {0,1,2}
            obj.(name) = []; 
          case 4
            obj.(name) = model.para(iNode-1).IAVpercent.(name)(1,resIdx); 
        end
      end
    end

%%
    function obj = Initialize(obj,model)
      numPara = model.numParaList(model.currentNode);
      rng(model.rngSetting); 

      randArray = rand(numPara.f4+numPara.f5+numPara.f6, model.numCase); 

      iRand = 1;
      for name = model.nameList.para
        switch obj.target.(name).flag
          case {0,1,2} 
            obj.(name) = []; %for save memory. Allocate later
          case 4
            obj.(name) = obj.ComputeParaFullRange(name,randArray(iRand,:),model.paraIAVpercent.def);
            iRand = iRand + 1;
          otherwise
        end 
      end
    end

%% 
    function y = allocateFlag0to2(obj,name,model)
      switch obj.target.(name).flag
        case {0,1}
          y = zeros(1,model.numCase);
        case 2
          y = model.paraIAVpercent.fix .* ones(1,model.numCase);
      end
    end

%%
    function y = ComputeParaFullRange(obj,name,randArray,paraRange)
      switch name
        case ''
        otherwise  
          y = paraRange .* randArray; %Find from uniform distribution [0, model.paraIAVpercent]
      end
    end

%% 
    function plotLine(obj,categoryList,statsNameList)
      figure;
      subPlotLayout = MyFunc.gridXYSubplot(categoryList);
      tscale = TspanScale(obj.tspan,[1750 2015]);
      for iCat = 1:length(categoryList)
        subplot(subPlotLayout(1),subPlotLayout(2),iCat);
        ylabel(MyFunc.nameLabel(categoryList(iCat),"Parameter")); hold on
        y = getErrorBar(obj.(categoryList(iCat)),statsNameList);

        if isempty(y)
          continue; 
        end

        errorbar(tscale.tspan,y(:,1),y(:,2),y(:,3),'-ko','MarkerFaceColor','auto'); 

        mod = ModifyPlotStyle;
        mod.isUseSameAxis = "off";
        mod.setXlimTimeSubplot([1750 2015],subPlotLayout(2));
      end
    end
  end
end

function y = getErrorBar(data,statsNameList)
  switch statsNameList(1)
    case 'ave', y = mean(data,2);
    case 'med', y = median(data,2);
  end

  prcXX = MyCalc.(statsNameList(2))(data);
  y(:,2) = abs(y(:,1)-prcXX(:,1));
  y(:,3) = abs(y(:,1)-prcXX(:,2));
end
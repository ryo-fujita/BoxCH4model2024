% Contains particle filter code for smoothing. 
% It extracts timeseries of posterior distributions of atmospheric components 
% or parameters (e.g., objects_SmoothSummaryPostFull_atmos.mat) 
% using their filtered distributions (e.g., objects_atmosSummary.mat).
% Ryo Fujita 2024

classdef Smoothing
  methods(Static)
    function objBest = loadPosterior(varName,model)
      filename = Smoothing.filenmPostFull(model,varName);

      if model.lagFilterTerm == model.numNode && exist(filename,'file')
        load(filename,'objBest'); 
      else
        error("Error: loadPosterior(varName,model)")
      end
    end

%% 
    function y = filenmPostFull(model,varName)
      y = strcat(model.workdir,'/objects_SmoothSummaryPostFull_',varName,'.mat');
    end

%% Run smoother using stored filtered data at each iNode 
    function objS = runFixedLagSmootherOffline(obj,varList,logicSaveMat)
    % argument information:
    % varList = "atmos", "para" or ["atmos","para"]

      if obj.numNode-1 < 2
        error('obj.numNode-1 < 2'); 
      end

      if ismember(obj.randType,["Prior","Posterior","PosteriorMean"])
        objS = []; 
        return; 
      end

      % load objects_*Summary.mat
      obj = obj.loadFilteredData(varList); 

      for iRun = 1:size(obj.judge,1)
        for iNodeS = 2:obj.numNode-1         
          % allocate stored original filtered data to objS
          if iRun == 1 && iNodeS == 2
            objS = obj;
          else
            for var = varList
              objS.(var)(iRun,iNodeS) = obj.(var)(iRun,iNodeS);
            end
          end

          if iNodeS <= obj.lagFilterTerm+1
            iNodeListS = iNodeS:-1:2;
          else
            iNodeListS = iNodeS:-1:iNodeS-obj.lagFilterTerm;
          end

          % reconstruct smoothed distributions at each iNode using original filtered data & index position info
          for iNode = iNodeListS
            objS = reconstructSmoothedData(obj,objS,iRun,iNode,iNodeS,varList);
          end
        end
        objS = allocateBlankToEdgeNode(obj,objS,iRun,iNode,varList);
      end

      if all(ismember(["atmos", "para"], varList)) && logicSaveMat
        atmos = objS.atmos;
        para  = objS.para;
        save(Smoothing.filenmSummary(obj),'atmos','para','-v7.3');

        Smoothing.extractPosteriorOffline("atmos",obj);
        Smoothing.extractPosteriorOffline("para",obj);
      end

      %
      function objS = reconstructSmoothedData(obj,objS,iRun,iNode,iNodeS,varList)
        if iNode ~= iNodeS
          for varName = varList
            objS.(varName)(iRun,iNode) = obj.judge(iRun,iNode).updateResults(...
              obj.judge(iRun,iNodeS).parentIdx,objS.(varName)(iRun,iNode)); 
          end
        end

        for varName = varList
          objS.(varName)(iRun,iNode) = obj.judge(iRun,iNode).updateResults(...
            obj.judge(iRun,iNodeS).resIdx,objS.(varName)(iRun,iNode)); 
        end
      end

      %
      function objS = allocateBlankToEdgeNode(obj,objS,iRun,iNode,varList)
        for varName = varList
          % iNode=1 or numNode is irrelevant to filtering, so leave blank data.
          objS.(varName)(iRun,1) = obj.judge(iNode).updateResults([],objS.(varName)(iRun,1));
          objS.(varName)(iRun,obj.numNode) = obj.judge(iNode).updateResults([],objS.(varName)(iRun,obj.numNode));
        end
      end
    end

%% 
    function objBest = extractPosteriorOffline(varName,model)
      obj = Smoothing.loadLagSmoothedData(model,varName).(varName);
      numRun    = size(obj,1);
      fieldList = model.nameList.(varName);
      
      switch varName
        case 'atmos', objBest = Atmosphere;
        case 'para',  objBest = ParameterSourceSink;
      end
      
      for iRun = 1:numRun
        objBest_each = Smoothing.extractEachRun(model,obj(iRun,:),fieldList);
        objBest = Smoothing.combineEachRun(...
          iRun,model,obj,fieldList,objBest_each,objBest);
      end

      if isempty(objBest.tspan), error('Empty smoothed posterior!'); end 
      save(Smoothing.filenmPostFull(model,varName),'objBest','-v7.3');
    end

%% 
    function objS = loadLagSmoothedData(model,varNameList)
      filename = Smoothing.filenmSummary(model);

      if exist(filename,'file') == 0
        objS = Smoothing.runFixedLagSmootherOffline(model,varNameList,false);
      elseif length(varNameList) == 1
        objS = load(filename,varNameList); 
      else
        objS = matfile(filename); 
      end   
    end

%%
    function y = filenmSummary(model)
      y = strcat(model.workdir,...
        '/objects_SmoothSummary_',Smoothing.lagName(model),'.mat');

      if model.issp ~= 0 && strcmp(model.randType,"OSSE")
        y = strrep(y,".mat",strcat(MyNameList.ssp(model.issp),"OSSE.mat"));
      elseif model.issp ~= 0
        y = strrep(y,".mat",strcat(MyNameList.ssp(model.issp),".mat"));
      end
    end

%% 
    function y = lagName(model)
      if model.lagFilterTerm == model.numNode
        y = "Full";
      else
        y = strcat('L',num2str(model.lagFilterTerm,'%02.0f'));
      end
    end

%%
    function objBest_each = extractEachRun(model,inputObj,fieldList)
      objBest_each.tspan = []; 
      for name = fieldList, objBest_each.(name) = [];
        if isempty(inputObj(model.numNode-1).(name))
          continue; 
        end

        for iNode = 1:size(inputObj,2)
          if iNode == 1 || iNode == model.numNode
            continue; 
          end 

          if strcmp(name,fieldList(1))
            objBest_each.tspan = Smoothing.combineObjectINode(...
              iNode,inputObj,inputObj(iNode).tspan.',objBest_each.tspan); 
          end 

          objBest_each.(name)  = Smoothing.combineObjectINode(...
            iNode,inputObj,inputObj(iNode).(name),objBest_each.(name));   
        end
      end
       
      if isa(inputObj,'Atmosphere') 
        idx = find(objBest_each.tspan>= 1970);
        tspan_new1970 = 1970:0.25:objBest_each.tspan(end);
        for name = fieldList
          obj_new1970 = interp1(objBest_each.tspan(idx), ...
            objBest_each.(name)(idx,:),tspan_new1970,'linear'); 
          objBest_each.(name)(idx(1):idx(1)+length(tspan_new1970)-1,:) = obj_new1970;

          if length(idx) > length(tspan_new1970)
            objBest_each.(name)(idx(1)+length(tspan_new1970):end,:) = [];
          end
        end

        objBest_each.tspan(idx(1):idx(1)+length(tspan_new1970)-1) = tspan_new1970;

        if length(idx) > length(tspan_new1970)
          objBest_each.tspan(idx(1)+length(tspan_new1970):end) = [];
        end    
      end  

      if isa(inputObj,'ParameterSourceSink') ...
        && isa(inputObj(1).IAVpercent,'ParaTimeVar') ...
        && ~isempty(inputObj(model.numNode-1).(name))  
        objBest_each = Smoothing.combineObjectINodeParaTimeVar(...
          inputObj,objBest_each,fieldList,model);
      end

      if isa(inputObj,'ParameterSourceSink') ...
        && isa(inputObj(1).hyper,'HyperParameter') ...
        && ~isempty(inputObj(model.numNode-1).(name))  
        objBest_each = Smoothing.combineObjectINodeHyper(...
          inputObj,objBest_each,fieldList,model); 
      end
    end

%% 
    function objBest = combineObjectINode(iNode,inputObj,objBestINode,objBest)
      iLast = size(objBestINode,1);

      switch class(inputObj) 
      case 'Atmosphere'
        if iNode == 2 
          objBest = objBestINode(1:iLast-1,:);
        elseif iNode == size(inputObj,2)-1
          objBest = [objBest; objBestINode];
        else
          objBest = [objBest; objBestINode(1:iLast-1,:)];
        end
      case 'ParameterSourceSink' 
        if iNode == 2 
          objBest = objBestINode(1:iLast-1,:);
        elseif iNode == size(inputObj,2)-1
          objBest = [objBest; objBestINode(2:iLast,:)];
        else
          objBest = [objBest;objBestINode(2:iLast-1,:)];
        end
      case {'ParaTimeVar','HyperParameter'} 
        objBest = [objBest;objBestINode];
      otherwise
        error(class(objBest) + " is not supported!")
      end
    end

%%
    function objBest_each = combineObjectINodeParaTimeVar(inputObj,objBest_each,fieldList,model)
      objBest_each.IAVpercent = ParaTimeVar;
      for iNode = 3:size(inputObj,2)-1 
        objBest_each.IAVpercent.tspan = Smoothing.combineObjectINode(...
          iNode,inputObj(1).IAVpercent,round(inputObj(iNode).tspan(1)),...
          objBest_each.IAVpercent.tspan); 
        for name = fieldList
          if isempty(inputObj(iNode).IAVpercent.(name))
            inputObj(iNode).IAVpercent.(name) = ...
              inputObj(iNode).IAVpercent.allocateFlag0to2(name,model); 
          end 

          objBest_each.IAVpercent.(name) = Smoothing.combineObjectINode(...
            iNode,inputObj(1).IAVpercent,inputObj(iNode).IAVpercent.(name),...
            objBest_each.IAVpercent.(name));
        end
      end  
    end

%%
    function objBest_each = combineObjectINodeHyper(inputObj,objBest_each,fieldList,model)
      objBest_each.hyper = HyperParameter;
      for iNode = 1:size(inputObj,2)-1 
        objBest_each.hyper.tspan = Smoothing.combineObjectINode(...
          iNode,inputObj(1).hyper,round(inputObj(iNode).tspan(1)),...
          objBest_each.hyper.tspan); 
        for name = model.nameParaList(iNode).hyper
          if inputObj(iNode).ismemberHyper(name)
            objBest_each.hyper.(name) = Smoothing.combineObjectINode(...
              iNode,inputObj(1).hyper,inputObj(iNode).hyper.(name),...
              objBest_each.hyper.(name));
          end
        end
      end  
    end

%%
    function inputObjBest_Comb = combineEachRun(...
      iRun,model,inputObj,fieldList,inputObjBest_each,inputObjBest_Comb)
      if iRun == 1
        inputObjBest_Comb.tspan = []; 
      end
 
      if isempty(inputObjBest_each) ...
        || size(inputObjBest_each.(fieldList(1)),1) ~= length(inputObjBest_each.tspan)
        return; 
      end

      if isempty(inputObjBest_each)
        return; 
      end
 
      if size(inputObjBest_each.(fieldList(1)),1) ~= length(inputObjBest_each.tspan)
        error(['The vector length of inputObjBest_each.(fieldList(1)) and ' ...
          'inputObjBest_each.tspan is not consistent!']); 
      end

      for name = fieldList
        inputObjBest_Comb.(name) = CombineObjectName(...
          iRun,model,name,fieldList,inputObj,inputObjBest_each,inputObjBest_Comb);
      end

      if isempty(inputObjBest_Comb.tspan)
        inputObjBest_Comb.tspan = inputObjBest_each.tspan; 
      end

      if isa(inputObj,'ParameterSourceSink') ...
        && isa(inputObj(1).IAVpercent,'ParaTimeVar') ...
        && ~isempty(inputObj(model.numNode-1).(name))
        inputObjBest_Comb.IAVpercent = ParaTimeVar;
        for name = fieldList
          inputObjBest_Comb.IAVpercent.(name) = CombineObjectName(...
            iRun,model,name,fieldList,inputObj(1).IAVpercent,...
            inputObjBest_each.IAVpercent,inputObjBest_Comb.IAVpercent);
        end
        if isempty(inputObjBest_Comb.IAVpercent.tspan)
          inputObjBest_Comb.IAVpercent.tspan = inputObjBest_each.IAVpercent.tspan; 
        end
      end

      if isa(inputObj,'ParameterSourceSink') ...
        && isa(inputObj(1).hyper,'HyperParameter') ...
        && ~isempty(inputObj(model.numNode-1).(name)) 
        inputObjBest_Comb.hyper = HyperParameter;
        for name = model.nameParaList(1).hyper
          inputObjBest_Comb.hyper.(name) = CombineObjectName(...
            iRun,model,name,model.nameParaList(1).hyper,...
            inputObj(1).hyper,inputObjBest_each.hyper,inputObjBest_Comb.hyper);
        end
        if isempty(inputObjBest_Comb.hyper.tspan)
          inputObjBest_Comb.hyper.tspan = inputObjBest_each.hyper.tspan; 
        end
      end  

      disp([''])
    end

%%
    function objS = fixedLagSmootherOnline(obj)
      if obj.currentNode == 1 || obj.currentNode == obj.numNode
        return; 
      end

      if obj.currentNode <= obj.lagFilterTerm+1  
        iNodeList = obj.currentNode:-1:2;
      else
        iNodeList = obj.currentNode:-1:obj.currentNode-obj.lagFilterTerm;
      end
 
      filename  = strcat('tmp_objSmooth_',num2str(obj.currentRun,'%05.0f'),'.mat');      
      if obj.currentNode == 2
        objS = obj;
      else
        objS = load(filename);
        objS.atmos(obj.currentNode) = obj.atmos(obj.currentNode);
        objS.para(obj.currentNode)  = obj.para(obj.currentNode);
      end

      for iNode = iNodeList
        if iNode ~= obj.currentNode
          objS.atmos(iNode) = obj.judge(iNode).updateResults(...
            obj.judge(obj.currentNode).parentIdx,objS.atmos(iNode));
          objS.para(iNode) = obj.judge(iNode).updateResults(...
            obj.judge(obj.currentNode).parentIdx,objS.para(iNode));
        end
        objS.atmos(iNode) = obj.judge(iNode).updateResults(...
          obj.judge(obj.currentNode).resIdx,objS.atmos(iNode));
        objS.para(iNode) = obj.judge(iNode).updateResults(...
          obj.judge(obj.currentNode).resIdx,objS.para(iNode));
      end
      atmos = objS.atmos;
      para  = objS.para;
      save(filename,'atmos','para','-v7.3'); 
    end

%%
    function y = checkValidLagTerm(lagTermList,model)
      i = 1;
      for lag = lagTermList
        if lag <= model.numNode
          y(i) = lag;
          i = i + 1;
        end
      end
    end

%%
    function lgdNameList = convertLagToNameList(lagTermList,model)
      lgdNameList = strings(1,length(lagTermList));
      for i = 1:length(lagTermList)
        model.lagFilterTerm = lagTermList(i);
        lgdNameList(i) = Smoothing.lagName(model);
      end
    end

%% 
    function numCase = countUniqueNumCase(lagList,varName,model)
      for lag = lagList
        model.lagFilterTerm = lag;
        lagName = Smoothing.lagName(model);
        objS    = Smoothing.loadLagSmoothedData(model,varName).(varName);
        numCase.(lagName) = struct('fbb',zeros(size(objS)),'Egeo',zeros(size(objS)));
        for iRun = 1:size(objS,1)
          for iNode = 1:size(objS,2)
            if isempty(objS(iRun,iNode).fbb)
              continue; 
            end
            numCase.(lagName).fbb(iRun,iNode)  = length(unique(objS(iRun,iNode).fbb(1,:)));
            numCase.(lagName).Egeo(iRun,iNode) = length(unique(objS(iRun,iNode).Egeo(1,:)));
          end
        end
      end
      numCase.yLimMax = model.numCase_org * size(objS,1);
    end

%%
    function y = convertStructToArray(numCase)
      lagNameList  = string(fieldnames(numCase)).';
      iNodeList    = 1:size(numCase.(lagNameList(1)).fbb,2);
      y = struct('fbb',zeros(length(lagNameList)-1,length(iNodeList)),...
        'Egeo',zeros(length(lagNameList)-1,length(iNodeList)));
      iLag = 0;
      for lagName = lagNameList
        if strcmp(lagName,"yLimMax"); continue; end
        iLag = iLag + 1;
        for iNode = iNodeList
          y.fbb(iLag,iNode)  = sum(numCase.(lagName).fbb(:,iNode));
          y.Egeo(iLag,iNode) = sum(numCase.(lagName).Egeo(:,iNode));
        end
      end
      y.yLimMax = numCase.yLimMax;
    end

%%
    function plotLagListVsNumCase(lagList,iNodeList,model,varargin)
      if ~isempty(varargin)
        subPlotLayout = varargin{1};
      else
        subPlotLayout = [3 4];
      end

      numCase     = Smoothing.countUniqueNumCase(lagList,"para",model);
      numCaseList = Smoothing.convertStructToArray(numCase);
      
      figure; iPlot = 0; iFig = 0;
      for iNode = iNodeList
        if iNode >= model.numNode
          break; 
        end
        iPlot = iPlot + 1;
        if iPlot > subPlotLayout(1)*subPlotLayout(2)
          iFig = iFig + 1;
          if iFig == 1
            modify = ModifyPlotStyle(model); 
          end
          modify.setSizeFig(subPlotLayout);
          filenm = strcat('NumCase_vs_Lag',num2str(lagList(1)),...
            '-',num2str(lagList(length(lagList))),'_',num2str(iFig));
          print(filenm,'-dpdf');
          figure; iPlot = 1; 
        end

        subplot(subPlotLayout(1),subPlotLayout(2),iPlot)
        if iNode == model.numNode -1
          title({strcat(num2str(round(model.paraTargetList(51).tspan(2)-1))),...
            strcat(" (iNode = ",num2str(iNode),")")}); 
          hold on
        else
          title({strcat(num2str(round(model.paraTargetList(iNode).tspan(1)))),...
            strcat(" (iNode = ",num2str(iNode),")")}); 
          hold on
        end

        plot(lagList,log10(numCaseList.fbb(:,iNode)),'-o','MarkerFaceColor','auto'); 
        hold on %flag4 
        plot(lagList,log10(numCaseList.Egeo(:,iNode)),'-o','MarkerFaceColor','auto'); %flag1 
        ylim([0 log10(numCaseList.yLimMax)]);

        if rem(iPlot,subPlotLayout(2)) == 1
          ylabel('Log10(Particle N)'); 
        end
        
        if floor(iPlot/subPlotLayout(2)) >= subPlotLayout(1)-1
          xlabel('Lag Node'); 
        end

        if iPlot == 1
          lgd = legend; lgd.String = fieldnames(numCaseList); 
        end
      end
      modify.setSizeFig(subPlotLayout);
      filenm = strcat('NumCase_vs_Lag',num2str(lagList(1)),...
        '-',num2str(lagList(length(lagList))),'_',num2str(iFig+1));
      print(filenm,'-dpdf');
    end

  end
end

%%
function out = CombineObjectName(...
  iRun,model,name,fieldList,inputObj,inputObjBest,inputObjBest_combined)
  if ~isa(inputObj,'ParameterSourceSink')
    if isempty(inputObjBest_combined.(name))
      out = inputObjBest.(name);
    else
      out = [inputObjBest_combined.(name),inputObjBest.(name)];
    end
    return
  end
  
  if iRun == 1
    if model.paraTargetList(length(model.paraTargetList)).(name).flag == 0
      out = inputObj(1,1).target.(name).def .* ones(size(inputObjBest.(fieldList(1))));
    else
      out = inputObjBest.(name);
    end
  else
    if model.paraTargetList(length(model.paraTargetList)).(name).flag == 0
      out = [inputObjBest_combined.(name), ...
        inputObj(1,1).target.(name).def .* ones(size(inputObjBest.(fieldList(1))))];
    else
      out = [inputObjBest_combined.(name), inputObjBest.(name)];
    end
  end
end

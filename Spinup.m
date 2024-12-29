% Handles the spin-up calculations.
% It keeps spinup under a condition at iNode = 2 in Test.m
% Ryo Fujita 2024

classdef Spinup
  properties
    type = "noise" % "uniform" or "noise" or "Cauchy" or []
    paraNoisePercent = 0.1 
    step = 0 
    needsCalcAtmosEquilibrium = true 
    exists   = false %initial: false & updated as true after spinup file has been created
    isSave   = true;
  end

  methods
    function obj = Spinup(timeSpinUpYr,randType)
      if timeSpinUpYr == 0 || ismember(randType,["Prior","Posterior","PosteriorMean"])
        obj.needsCalcAtmosEquilibrium = false; 
      end 
    end

%% 
    function obj = saveMatFileTmpParaTarget(obj,model)
      iNode = model.currentNode;
      paraTarget = model.paraTargetList(1);

      switch obj.type
        case "uniform"
          for name = model.nameList.para
            paraRange = paraTarget.(name).max - paraTarget.(name).min; 
            paraTarget.(name).min = ...
              max(min(model.para(iNode).(name)(end,:))-paraRange/50, paraTarget.(name).min); 
            paraTarget.(name).max = ...
              min(max(model.para(iNode).(name)(end,:))+paraRange/50, paraTarget.(name).max);
          end
          obj.needsCalcAtmosEquilibrium = false; 
        case {"noise","Cauchy"} 
          for name = model.nameList.para
            paraTarget.(name).def = model.para(2).(name)(end,:);
          end
        otherwise
          error("Not yet developed")
      end

      save(strcat(model.workdir,'/tmp_paraTarget_Spinup_initial_',...
        num2str(model.currentRun,'%05.0f'),'.mat'),'paraTarget');
      obj.exists = true;
    end

  end

%%
  methods(Static)
%%
    function saveMatFileTmpPosterior(model)
      filename = model.tmpMatFileName(model.currentRun,'tmp_objectsSpinUp_');
      model = model.ClearMemorySelect;
      save(filename,'model');
    end

%% 
    function CreateSummaryMat(model)
      for varName = ["atmos","para","judge"]
        for iRun = 1:model.numRun
          filename  = model.tmpMatFileName(iRun,'tmp_objectsSpinUp_');

          if iRun == 1
            var = load(filename).model.(varName)(2); %extract iNode = 2 only
            var(model.numRun) = var;
          else
            var(iRun) = load(filename).model.(varName)(2);  
          end  

          disp(strcat("Finish reading ",filename, " for ", varName," ",datestr(datetime('now'))))
        end
        Info = ReadModelParameter; 
        save(strcat('objectsSpinUp_',varName,'Summary.mat'),'Info','var','-v7.3'); 
        clear var
      end
      disp(strcat("obj.CreateSpinUpSummaryMat Finish ",datestr(datetime('now'))))
    end

%% 
    function model = loadMatFile(model)
      for propName = ["judge","atmos","para"]
        filename = strcat("objects_",propName,"SpunUpSummary.mat");
        model = model.loadMatFile(propName,filename);
      end
      model.currentIterNum = model.currentIterNum + 1; 
    end

  end
end
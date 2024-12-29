% The primary program responsible for controlling 
% model simulation, smoothing processes, and generating outputs. 
% It can be run as RunModel().main. 
% Ryo Fujita 2024
classdef RunModel
  properties 
    setting = "normal"  
    varNameList = ["atmos","para","judge"]
    numBatch = 1 %default: 1 in "normal" 
    numLoop
    waitTime = 60 %default: 60 seconds
    PBS_JOBNAME
    workdir = "."
    workdirList 
  end

  properties (Dependent)
    numRun
    matfileList 
  end

  methods
    function obj = RunModel(varargin)
      if ~isempty(varargin) 
      for i = 1:length(varargin)
        if isempty(varargin{i}), continue; end
        switch i
          case 1, obj.setting = varargin{1};
          case 2, obj.PBS_JOBNAME = varargin{2};
          case 3, obj.workdir = varargin{3};
        end
      end
      end

      info            = ReadModelParameter(obj.workdir); 
      obj.numLoop     = info.numLoop;
      obj.workdirList = info.workdirList;

      if ~isempty(info.numBatch)
        obj.numBatch = info.numBatch;
      end
    end

%% Executed as RunModel(PBS_JOBNAME).main 
    function main(obj)
      % Clean up current directory
      fclose('all'); close all; 
      cleanupDirectory();

      % Set variables
      addpath constants emission 
      model = Model(obj.workdir);
      printInpfileInfo(model);

      % Run BoxModel
      obj.RunBoxModel;

      % Combine ensemble simulations
      obj.CreateSummaryMat(model);

      % Smoothing
      model.currentNode = model.numNode; % need to set for Smoothing 
      model.lagFilterTerm = model.numNode; % need to set for Smoothing
      Smoothing.runFixedLagSmootherOffline(model,["atmos","para"],true); 

      % Write filtered results
      %model = model.loadFilteredData(obj.varNameList);
      %writeFilteredResults(model,1:model.numNode);

      % Write summary posterior stats
      runWriteStatsPosterior(obj.workdir); 
      clear model 

      % Create main figures using multiple inventory scenario means
      if isReadyMultiScenario(obj)
        model = Model(obj.workdirList(1),"Posterior",9);
        Fig1(obj.workdirList,model); close all;
        Fig6(obj.workdirList,model); close all;
    
        runPlot = RunPlotBest(obj.workdirList);
        runPlot.ParaAllSuppliment; close all; 
        runPlot.Prctile; close all; 
        clear runPlot

        runWriteStatsPosterior(obj.workdirList);
      end  

      %Create figures for each prior scenario (optional)
      runPlot = RunPlotBest(obj.workdir); 
      runPlot.Prctile;   close all;  
      runPlot.ParaTimeVar; close all; 
      runPlot.ParaHyper;  
      close all; 
      clear runPlot

      disp('RunModel Finish')
      delete tmp*.mat;
    end

%% 
    function y = isReadyMultiScenario(obj)
      y = isPWDLastDir(obj) && isReadyOtherDir(obj);
    end

%% 
    function y = isPWDLastDir(obj)
      if isempty(obj.workdirList) || length(obj.workdirList) == 1
        y = false; 
        disp("False: isPWDLastDir(obj)")
        disp("obj.workdirList: " + obj.workdirList)
        disp("length(obj.workdirList): " + num2str(length(obj.workdirList)))
        disp("")
        return; 
      end
 
      y = strcmp(obj.workdir,obj.workdirList(end));
    end

%% 
    function y = isReadyOtherDir(obj)
      for iDir = 1:length(obj.workdirList) 
        workDir = obj.workdirList(iDir);
        count = 0;
        while ~exist(strcat(workDir,'/objects_atmosSummary.mat'),'file') ...
              || ~exist(strcat(workDir,'/objects_paraSummary.mat'),'file') ...
              || ~exist(strcat(workDir,'/objects_judgeSummary.mat'),'file') 
          pause(60)
          count = count + 1;

          disp(strcat([datestr(datetime('now')),' paused 60s until ', ...
            strcat(workDir,'/objects_atmosSummary.mat'),' & ', ...
            strcat(workDir,'/objects_paraSummary.mat'),' is created']))

          assert(count < 180,strcat([datestr(datetime('now')), ...
            strcat(workDir,'/objects_atmosSummary.mat'),' & ', ...
            strcat(workDir,'/objects_paraSummary.mat'), ...
            ' do not exist over the 3 hr. Please check if there are some errors']))
        end
      end
      y = true;
    end

%% 
    function RunBoxModel(obj)  
      switch obj.setting
      case 'normal'  
        for iRun = 1:obj.numRun
          disp(strcat(['iRun = ',num2str(iRun,'%05d')]))
          Model(obj.workdir).main(iRun); 
          clear ans 
        end

      % Contact Ryo Fujita if you would like to use the batch process code
      case 'batch' 
        c = parcluster('local'); disp(c)
        ctrJob = ControlBatchJob(obj.PBS_JOBNAME,obj.workdir); 

        for iLoop=1:obj.numLoop
          disp(datetime('now')); 
          findJob(c) 
          ctrJob = ctrJob.WaitJob();  
          job    = ctrJob.SubmitJob(c,iLoop,obj.numBatch);  
          pause(obj.waitTime); 
          ctrJob = ctrJob.DeleteJob(job,iLoop,obj.numBatch); 
        end
      end
      disp(strcat("obj.RunBoxModel() Finish ",datestr(datetime('now'))))
    end

%% 
    function CreateSummaryMat(obj,model)
      for iVar = 1:length(obj.varNameList)        
        for iRun = 1:obj.numRun
          varName   = obj.varNameList(iVar);
          filename  = model.tmpMatFileName(iRun);

          if ~exist(filename,'file'), obj.RunMainAgain(iRun,filename); end  

          if iRun == 1
            var = load(filename).obj.(varName);
            var(obj.numRun,:) = var;
          else
            var(iRun,:) = load(filename).obj.(varName);  
          end  

          disp(strcat("Finish reading ",filename, " for ", varName," ",datestr(datetime('now'))))
        end

        % Save object into matfile 
        Info = ReadModelParameter; 
        outMatfile = strcat(model.workdir,"/",obj.matfileList(iVar));
        
        % if model.issp ~=0 && strcmp(model.randType,"OSSE")
        %   outMatfile = strrep(outMatfile,".mat",strcat(MyNameList.ssp(model.issp),"OSSE.mat"));
        % elseif model.issp ~=0
        %   outMatfile = strrep(outMatfile,".mat",strcat(MyNameList.ssp(model.issp),".mat"));
        % end

        save(outMatfile,'Info','var','-v7.3');
        clear var
      end
      disp(strcat("obj.CreateSummaryMat Finish ",datestr(datetime('now'))))
    end

%% 
    function RunMainAgain(obj,iRun,filename)
      disp(strcat([filename,' does not exist! Run Model(),main again.']))
      disp(strcat(['iRun = ',num2str(iRun,'%05d')]))
      Model(obj.workdir).main(iRun); 
    end

%%
    function numRun = get.numRun(obj)
      if strcmp(obj.setting,"normal")
        numRun = obj.numLoop; 
      end

      if strcmp(obj.setting,"batch")
        numRun = obj.numLoop * obj.numBatch; 
      end
    end

%%
    function matfileList = get.matfileList(obj)
      matfileList = strings(1,length(obj.varNameList));
      for iVar = 1:length(obj.varNameList)
        varName = obj.varNameList(iVar);
        matfileList(iVar) = strcat('objects_',varName,'Summary.mat');
      end
    end

  end
end      

%% 
function cleanupDirectory()
  delete objects_*Summary.mat tmp*.mat tmp*.txt tmp*.png tmp*.fig;
  delete judge_*.txt
end

%% 
function printInpfileInfo(model)
  fprintf("%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n",...
    strcat(['Main Input Files are as below: ',...
    model.inputFile,model.emsCH4File])); 
end


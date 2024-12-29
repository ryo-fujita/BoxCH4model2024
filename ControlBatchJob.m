classdef ControlBatchJob
  properties
    infile = '../../jobStatus.txt'
    status
    jobName
    workdir 
    jobNameRunning 
    jobNameNextRun
    qsubTime
    qsubHistory
  end
  
  methods
    function obj = ControlBatchJob(jobName,workdir)
      if nargin > 0
        obj.jobName = jobName;
        obj.workdir = workdir;
        if strcmp(extractBefore(obj.jobName,2),"1")
          fileID = fopen(obj.infile,'w');
          fprintf(fileID,'%s\t%s\t%s\n',obj.jobName,datestr(datetime('now')),'has started'); 
          obj.status  = 'ready';
          fclose(fileID);
          obj = obj.WriteJobStatus('ready');
        else
          obj.status  = 'wait';
          pause(60); 
        end 
        obj.qsubHistory = obj.setQsubHistory; 
        obj.qsubTime    = obj.setQsubTime; 
      end
    end

%%
    function obj = WriteJobStatus(obj,status)
      fileID = fopen(obj.infile,'a');
      obj.status = status;
      switch obj.status
        case 'ready'
          fprintf(fileID,'%s\t%s\t%s\n',obj.jobName,datestr(datetime('now')),'*Ready to submit jobs');
          disp(strcat([datestr(datetime('now')),' now c.Jobs is empty and ready to submit']))
        case 'submit' 
          fprintf(fileID,'%s\t%s\t%s\n',obj.jobName,datestr(datetime('now')),'*Submit jobs (start)');
        case 'delete'
          fprintf(fileID,'%s\t%s\t%s\n',obj.jobName,datestr(datetime('now')),'*Delete jobs (start)');
      end
      
      fclose(fileID);
    end
    
%% 
    function y = readEndLine(obj,type)
      fid = fopen(obj.infile);
      while ~feof(fid)
        tline = fgetl(fid); 
      end
    
      if tline == -1
        y = '';   
        fclose(fid); 
        return; 
      end

      switch type
        case 'jobName'
          lines = split(tline);
          y = lines{1};
        case 'status'
          y = extractAfter(tline,'*');
        case 'all'
          y = tline;
      end
      fclose(fid);
    end
    
%% 
    function  obj = WaitJob(obj)
      count = 0;
      while ~strcmp(obj.readEndLine('status'),'Ready to submit jobs') || ...
            ~strcmp(obj.jobName,obj.jobNameNextRun)
        pause(60)
        disp(strcat([datestr(datetime('now')),' paused 60s until jobStatus shows ready and obj.jobNameNextRun is',obj.jobName])) %disp(c.Jobs)
        count = count + 1;
        assert(count < 720,strcat([datestr(datetime('now')),' WaitJob for 12 hr. SC System would be very busy or there would be some errors.']))
      end

      obj = obj.WriteJobStatus('submit');
    end

%%
    function y = get.jobNameRunning(obj)
      y = obj.readEndLine('jobName');
    end

%%
    function y = get.jobNameNextRun(obj)
      if strcmp(obj.jobName,obj.jobNameRunning)
        disp('This jobName is currently running!'); 
        y = obj.jobName; 
        return; 
      end

      obj.qsubHistory = obj.setQsubHistory; 

      for i = 1:length(obj.qsubHistory.jobName)
        qsubName = obj.qsubHistory.jobName{i};

        if strcmp(qsubName,obj.jobNameRunning)
          y = obj.qsubHistory.jobName{i+1}; 
          break
        else 
          y = obj.jobName; 
        end
      end
    end

%%
    function obj = DeleteJob(obj,job,iLoop,numBatch)
      timeout = 60*60*2; %2hr
      for iBatch=1:numBatch
        iRun=numBatch*(iLoop-1)+iBatch;
        jobNum = strcat("i",num2str(iRun,"%05d"));        
        ok = wait(job.(jobNum),'finished',timeout);
        if ok == true
          disp(strcat(jobNum,[' has been finished at ', datestr(datetime('now')), ', Job ID = ',num2str(job.(jobNum).ID)]))
          disp(job.(jobNum))   
          obj.findError(job,iRun); 
          if iBatch == 1, obj = obj.WriteJobStatus('delete'); end
          delete(job.(jobNum)); 
          job = rmfield(job,jobNum);
        else
          disp(strcat(jobNum,[' seems not be working properly. Skip to the next job.', datestr(datetime('now'))]))
        end
      end

      obj = obj.WriteJobStatus('ready');
      pause(61) 
    end

%%
    function y = isExcecuted(obj)
      y = exist(obj.infile,'file') ~= 0;
    end  

%%
    function y = setQsubTime(obj)
      numLine = length(obj.qsubHistory.jobName);
      if strcmp(obj.jobName,obj.qsubHistory.jobName(numLine))
        y = obj.qsubHistory.time(numLine);
      else
        error('Inconsistent jobName in file qsubHistory.txt & PBS_JOBNAME.txt');
      end
    end

%%
    function job = SubmitJob(obj,objCluster,iLoop,numBatch)
      for iBatch=1:numBatch
        iRun = numBatch*(iLoop-1)+iBatch;
        jobNum = strcat("i",num2str(iRun,"%05d"));

        job.(jobNum) = batch(objCluster,@RunModelForBatch,0,{obj.workdir,iRun},'CaptureDiary',true); pause(1); %2020.6.17
        disp(strcat(jobNum, [' has been submitted at ', datestr(datetime('now')), ', Job ID = ',num2str(job.(jobNum).ID)]))
      end
    end
  end
  
  methods(Static)
     function y = setQsubHistory()
      fid1 = fopen('./qsubHistory.txt','r');
      InputText = textscan(fid1, '%s %n','delimiter', ' ');
      fclose(fid1);  

      y.jobName  = InputText{1,1};
      y.time     = InputText{1,2};
     end

%%
    function findError(job,iRun)
      jobNum = strcat("i",num2str(iRun,"%05d"));
      if strcmp(job.(jobNum).State,'unavailable') 
        errfile = fopen('./error.log','a');
        fprintf(errfile,"%s %s\n",strcat(jobNum,[datestr(datetime('now')),' Job ID = ',num2str(job.(jobNum).ID),' unavailable']));
        fclose(errfile);

        errfile = fopen('./errorJobList.txt','a');
        fprintf(errfile,"%d\t%s\n",iRun,strcat(datestr(datetime('now')),' Job ID = ',num2str(job.(jobNum).ID),' unavailable'));
        fclose(errfile);            
      elseif ~isempty(job.(jobNum).Tasks(1).Error)
        errfile = fopen('./error.log','a');
        fprintf(errfile,"%s",getReport(job.(jobNum).Tasks(1).Error));
        fclose(errfile);

        errfile = fopen('./errorJobList.txt','a');
        fprintf(errfile,"%d\t%s\n",iRun,strcat(datestr(datetime('now')),' Job ID = ',num2str(job.(jobNum).ID),' Error'));
        fclose(errfile);
      else
      end
    end  
  end  
  
end


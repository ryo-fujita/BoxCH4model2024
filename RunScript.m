% Can be called by external script (e.g. shell script)
% when executing the model as batch processing. 
% RunModel.m can then be executed through this script.
% Ryo Fujita 2024
% 
% Arguments:
% 1. workdir: full path directory name running the model
% 2. PBS_JOBNAME: arbitrary job name when submitting this script. This will be used in ControlBatchJob.m 
% For more information to utilize the batch mode, please contact the author.

function RunScript(workdir,PBS_JOBNAME)
  workdir = string(workdir); 
  DateBegin = datetime('now'); disp(DateBegin)
  close all
  
  RunModel("batch",PBS_JOBNAME,workdir).main(); 
 %RunModel("normal",PBS_JOBNAME,workdir).main(); 

  disp('RunScript Finish')
  DateEnd = datetime('now'); disp(DateEnd)

  system(sprintf('echo "start:\t"%s >  %s/finish.txt\n', DateBegin,workdir)); 
  system(sprintf('echo "finish:\t"%s >> %s/finish.txt\n',DateEnd,  workdir));
  quit
end
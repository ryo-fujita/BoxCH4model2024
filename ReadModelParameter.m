% Handles the loading of basic properties from ModelInputParameter.info.
% You can add any properties of Model.m here, which you want to read through ModelInputParameter.info
% Ryo Fujita 2024

classdef ReadModelParameter 
   properties
    workdir
    inpfile = 'ModelInputParameter.info'
    numCase
    numLoop
    numBatch

    emsNatr = '' 
    emsBB   = '' 
    lossOH  = '' 

    inventory
    atmosTargetFile

    workdirList
    dirFile = 'dirname.txt'
  end

  methods
%%
    function obj = ReadModelParameter(varargin)
      if nargin == 0
        obj.workdir = '.';
      else
        obj.workdir = varargin{1};
      end

      obj = obj.readInfo;
      obj = obj.readDirList;
    end

%%
    function obj = readInfo(obj)
      [fid,errmsg] = fopen(strcat(obj.workdir,'/',obj.inpfile),'r'); 
      if ~isempty(errmsg); error(strcat(errmsg,strcat(obj.workdir,'/',obj.inpfile))); end

      while ~feof(fid)
        rowCell = textscan(fid,'%s',1,'delimiter','\n');
        nameRow = extractBefore(rowCell{1}{1},":");

        if ismember(nameRow,string(fieldnames(obj)))
          obj.(nameRow) = strtrim(extractAfter(rowCell{1}{1},":"));
        else
          disp(strcat(nameRow," is not properties. Please check the name is valid."))
          continue 
        end

        switch nameRow
          case {'numCase','numLoop','numBatch'}
            obj.(nameRow) = str2double(obj.(nameRow));
        end
      end
      fclose(fid);
    end

%%
    function obj = readDirList(obj)
      filename = strcat(obj.workdir,'/',obj.dirFile);
      pause(0.1) 
      if exist(filename,'file')==0
        obj.workdirList = obj.workdir; 
        return; 
      end
 
      if strcmp(version,'9.5.0.944444 (R2018b)')
        [fid,errmsg] = fopen(filename);
        if ~isempty(errmsg)
          error(strcat(errmsg,": ",filename)); 
        end

        C = textscan(fid,'%s\n');
        obj.workdirList = strings(length(C{1}),1);
        for i = 1:length(C{1})
          obj.workdirList(i) = C{1}{i};
        end
        fclose(fid);
      else
        obj.workdirList = readmatrix(...
          filename,'FileType','text','OutputType','string','Delimiter','\n');
      end

      if size(obj.workdirList,1) ~= 1 && size(obj.workdirList,1) > 1 
        obj.workdirList = obj.workdirList.'; %unify as column vector
      end
    end
  end
end
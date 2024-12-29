%
classdef ReadTargetFile
  properties
    workdir
    inputFile
    yrInitial
    timeSpinUpYr
    atmosOffset

  end

  methods
    function obj = ReadTargetFile(workdir,inputFile,yrInitial,timeSpinUpYr,atmosOffset)
      if nargin == 5
        obj.workdir = workdir;
        obj.inputFile = inputFile;
        obj.yrInitial = yrInitial;
        obj.timeSpinUpYr = timeSpinUpYr;
        obj.atmosOffset = atmosOffset;
      else
        error("Please set the correct nargin!")
      end
    end

%%
    function [nameList, atmosTargetList, paraTargetList,...
              numNode, useDistPrior, useHyperParameter,...
              namePara,numPara] = main(obj)

      % Read input file
      S = readlines(obj.inputFile,'EmptyLineRule','skip');
      rowHyphen = find(contains(S,lineBoundary + '-'));
      columnInfoStr = S(rowHyphen(1)+1:rowHyphen(2)-1);
      dataInfoStr   = S(rowHyphen(3)+1:end);

      % Get information 
      [term,atmosInfo,paraInfo] = getColumnInfo(columnInfoStr,obj.inputFile);
      
      if term.num ~= length(dataInfoStr)-2
        error("Please set correct term number"); 
      end

      if any(ismember(fieldnames(paraInfo).','distPrior'))
        useDistPrior = true; 
      else
        useDistPrior = false;
      end

      if any(ismember(fieldnames(paraInfo).','setHyper')) ...
        && any(cell2mat(struct2cell(paraInfo.setHyper)))
        useHyperParameter = true;
      else
        useHyperParameter = false;
      end

      nameList = getNameList(dataInfoStr,term,atmosInfo,paraInfo);
      
      if ~isConsistent(nameList,term,atmosInfo,paraInfo)
        error('The header information is invalid'); 
      end

      [atmosTargetList, paraTargetList] = obj.getDataInfo(dataInfoStr,term,atmosInfo,paraInfo);
      
      numNode  = size(atmosTargetList,2) + 2; 

      [namePara,numPara] = extractNameNumParaAllTerm(...
        nameList,paraTargetList,numNode,useHyperParameter);
    end

%% 
    function [atmosTarget, paraTarget] = getDataInfo(obj,inpStr,term,atmos,para)
      istColAtmos = 3;         
      istColPara  = 3+sum(cell2mat(struct2cell(atmos.numCol)));
      dataStr     = split(inpStr(istColAtmos:end),char(9));

      if size(dataStr,2) == 1
        dataStr = dataStr.'; 
      end
  
      yrTarget    = obj.convertToNum(dataStr(:,2)).';
      fieldname   = getFieldName(inpStr,atmos,para); 

      atmosTarget = obj.convertToTargetArray(atmos,dataStr,term,yrTarget,fieldname,istColAtmos); 
      paraTarget  = obj.convertToTargetArray(para,dataStr,term,yrTarget,fieldname,istColPara);
    end

%% 
    function targetArray = convertToTargetArray(...
      obj,tarInfo,dataStr,term,yrTarget,fieldname,ist)

      for nameID = tarInfo.name
        dataStrTmp = dataStr(:,ist:ist+tarInfo.numCol.(nameID)-1);
        tarData.(nameID) = zeros(size(dataStrTmp));
        for iCol = 1:size(dataStrTmp,2)
          tarData.(nameID)(:,iCol) = str2double(dataStrTmp(:,iCol));
          
          % Target files that contain special characters ($, @, etc.) and could not be 
          % converted directly to numerical values are read in individually
          if ~all(~isnan(tarData.(nameID)(:,iCol))) 

            % line numbers containing special characters
            for idx = find(isnan(tarData.(nameID)(:,iCol))).' 
              tarData.(nameID)(idx,iCol) = obj.convertToNum(dataStrTmp(idx,iCol),yrTarget(idx),nameID,tarData.(nameID)(idx,:));
            end
          end
        end

        ist = ist+tarInfo.numCol.(nameID);
      end

      yrTarget  = [obj.yrInitial - obj.timeSpinUpYr;yrTarget]; 
      targetArray(1,term.num) = struct('termName',"",'tspan',[]);
      
      for iTerm = 1:term.num
        targetArray(iTerm).termName = term.name(iTerm);
        switch tarInfo.type
          case "atmos", targetArray(iTerm).tspan = [yrTarget(iTerm),     yrTarget(iTerm+1)];
          case "para",  targetArray(iTerm).tspan = [yrTarget(iTerm)-0.5, yrTarget(iTerm+1)+0.5];
        end
        
        for name = string(fieldnames(tarData)).'
          for iField = 1:length(fieldname.(name))
            fname = fieldname.(name)(iField);
            targetArray(iTerm).(name).(fname) = tarData.(name)(iTerm,iField);
          end
          
          if strcmp(tarInfo.type,"para")  
            for fname = string(fieldnames(tarInfo).')
              switch fname
                case {'distPrior','setHyper'}
                  targetArray(iTerm).(name).(fname) = tarInfo.(fname).(name);
              end
            end
          end
        end
      end

      if strcmp(tarInfo.type,"atmos")
        targetArray = setTargetMinMax(targetArray);
      end

      %
      function atmosTargetList = setTargetMinMax(atmosTargetList)
        for i = 1:length(atmosTargetList)
        for gas = ["CH4","d13C","dD","D14C"]
          if ~ismember("min",fieldnames(atmosTargetList(i).(gas))) || isnan(atmosTargetList(i).(gas).min)
            atmosTargetList(i).(gas).min = atmosTargetList(i).(gas).ave - 3*atmosTargetList(i).(gas).sdev;
          end
          if ~ismember("max",fieldnames(atmosTargetList(i).(gas))) || isnan(atmosTargetList(i).(gas).max)
            atmosTargetList(i).(gas).max = atmosTargetList(i).(gas).ave + 3*atmosTargetList(i).(gas).sdev;          
          end
        end
        end
      end

    end

%%
    function y = convertToNum(obj,dataStrList,varargin)
      y = zeros(1,length(dataStrList));
      for i = 1:length(dataStrList)
        y(i) = obj.getVal(dataStrList(i),varargin{:});
      end
    end

%% 
    function y = getVal(obj,dataStrOrg,varargin)
      pat = ("@"+lettersPattern+"."+lettersPattern) | "@"+lettersPattern | "$"+lettersPattern; 
      dataStrOrgList = extract(dataStrOrg,pat);
      valOrgList = convertStr2Val(dataStrOrgList);
      y = eval(replace(dataStrOrg,dataStrOrgList,string(valOrgList)));

      function valList = convertStr2Val(dataStrList)
        valList = zeros(size(dataStrList));
        for i = 1:length(dataStrList)  
          dataStr = dataStrList(i);
          switch extract(dataStr,1) 
            case '$' %variable
              valList(i) = eval(replace(dataStr,"$","obj."));  
            % case '@' %function or expression  *only available to atmos
            %   yrTarget = varargin{1};
            %   gasID  = varargin{2};
            %   if size(yrTarget,1) ~= 1, error('Only scalar yr is available!'); end
            %   valList(i) = convertTarget(extractAfter(dataStr,1),yrTarget,gasID); 
          end
        end
      end

      
      % function y = convertTarget(str,yrTarget,gasID)
      %   switch str
      %     case 'Target' 
      %       field = "ave";
      %     case {'Target.min','Target.max','Target.minOffset','Target.maxOffset'}
      %       field = extractAfter(str,'.');
      %     otherwise
      %       field = str;
      %   end
      %   y = obj.ComputeAtmosTarget(gasID,field,yrTarget);
      % end

    end

%%
%     function val = ComputeAtmosTarget(obj,gas,field,tspan)
%       if ismember("Glb",fieldnames(obj.obs.(gas))) || contains(field,"OSSE")
%         t1 = 1750 <= tspan   & tspan <=  2015;
%         i1 = find(t1);
%         if isempty(i1), error('Invalid input tspan!'); end
% 
%         switch field
%           case {'ave','sdev','min','max'}
%             val(i1) = interp1(obj.obs.tspan.(gas).Glb,obj.obs.(gas).Glb.(field),tspan(i1),'linear');
%           case {'aveOSSE','sdevOSSE'}
%             assert(~isempty(obj.obsOSSE))
%             field = extractBefore(field,'OSSE');
%             val(i1) = interp1(obj.obsOSSE.tspan,obj.obsOSSE.(gas).(field),tspan(i1),'linear');
%           otherwise
%             error(strcat("There is no such a field in ",gas," target!"))
%         end
% 
%         if any(obj.atmosOffset.(gas)) && ~strcmp(field,'sdev')
%           val = val + obj.atmosOffset.(gas); 
%         end 
% 
%       else
%         switch gas
%           case "CH4", val = obj.ComputeAtmosTargetCH4(field,tspan); 
%           otherwise, error("Not prepared!"); 
%         end
%       end
% 
%       if ~all(isfinite(val),'all'), error(strcat([num2str(tspan),': val contains NaN!'])); end
%     end   
% 
% %%
%     function [val,growth] = ComputeAtmosTargetCH4(obj,field,tspan)
%       t1 = 0 < tspan & tspan <  1984.5;
%       t2 = 1984.5   <= tspan & tspan <= 2020.5;
%       i1 = find(t1);
%       i2 = find(t2);
% 
%       switch field
%         case 'ave'          
%           if sum(i1) >= 1
%             val(i1)    = interp1(obj.obs.tspan.CH4.M17, obj.obs.CH4.M17.global, tspan(t1),'linear');
%             growth(i1) = interp1(obj.obs.tspan.CH4.M17, gradient(obj.obs.CH4.M17.global), tspan(t1),'linear');
%           end
%           if sum(i2) >= 1
%             val(i2)    = interp1(obj.obs.tspan.CH4.NOAA,obj.obs.CH4.NOAA.global,tspan(t2),'linear');
%             growth(i2) = interp1(obj.obs.tspan.CH4.NOAA,gradient(obj.obs.CH4.NOAA.global),tspan(t2),'linear');
%           end
%         otherwise
%           error('There is no such a field in CH4 target!')
%       end
% 
%       if isempty(i1) && isempty(i2) 
%         error('Invalid input tspanEms!')
%       end
% 
%       if ~all(isfinite(val),'all'), error('val contains NaN!'); end
%     end   

  end
end

%%
function nameList = getNameList(inpStr,term,atmos,para)
  nameList.atmosConc = ["CH4","C13","CHD","C14"];
  for i = 1:term.num,  nameList.term(i) = string(textscan(inpStr(2+i),'%s',1)); end 

  tmpList        = strsplit(strtrim(inpStr(1)),{'\t'}); 
  nameList.atmos = tmpList(3:atmos.num+2);
  nameList.para  = tmpList(atmos.num+3:atmos.num+para.num+2);
end

%%
function fieldname = getFieldName(inpStr,atmos,para)
  tmpList = strsplit(strtrim(inpStr(2)),{'\t'});
  ist = 1;
  for i = 1:atmos.num
    fieldname.(atmos.name(i)) = tmpList(ist:ist+atmos.numCol.(atmos.name(i))-1); 
    ist = ist+atmos.numCol.(atmos.name(i));
  end
  for i = 1:para.num
    fieldname.(para.name(i)) = tmpList(ist:ist+para.numCol.(para.name(i))-1); 
    ist = ist+para.numCol.(para.name(i));
  end
end

%%
function y = isConsistent(nameList,term,atmos,para)
  y = isequal(nameList.term,term.name) & isequal(nameList.atmos,atmos.name) & isequal(nameList.para,para.name);
end

%% Read header info
function [term,atmos,para] = getColumnInfo(inpStr,inpfile)
  term.type = "term"; 
  term.num  = str2double(extract(inpStr(1),digitsPattern));                %% Read Term header information
  term.name = strsplit(strtrim(inpStr(2)),{'\t'}); 

  atmos.type = "atmos";
  atmos.num  = str2double(extract(inpStr(3),digitsPattern));               %% Read Atmos header information
  atmos.name = strsplit(strtrim(inpStr(4)),{'\t'});
  numCol     = str2double(strsplit(strtrim(inpStr(5)),{'\t'}));
  for i = 1:atmos.num, atmos.numCol.(atmos.name(i)) = numCol(i); end

  para.type = "para";
  para.num  = str2double(extract(inpStr(6),digitsPattern));                %% Read Para header information
  para.name = strsplit(strtrim(inpStr(7)),{'\t'}); 
  numCol    = str2double(strsplit(strtrim(inpStr(8)),{'\t'}));
  for i = 1:para.num, para.numCol.(para.name(i)) = numCol(i); end
  
  if contains(inpfile,'Base8')
    distPrior =  strsplit(strtrim(inpStr(9)),{'\t'});
    setHyper  =  str2double(strsplit(strtrim(inpStr(10)),{'\t'}));
    for i = 1:para.num
      para.distPrior.(para.name(i)) = distPrior(i);
      para.setHyper.(para.name(i))  = setHyper(i);
    end
  end
end

%%
function  [namePara,numPara] = extractNameNumParaAllTerm(nameList,paraTargetList,numNode,useHyperParameter)    
  namePara = struct('termName',"",'all',"",'f0',"",'f1',"",'f2',"",'f3',"",'f4',"",'f5',"",'f6',"",'f7',"",'f8',"",'f9',"",'fDynamic',"",'hyper',"");
  numPara  = struct('termName',"",'all',0,'f0',0,'f1',0,'f2',0,'f3',0,'f4',0,'f5',0,'f6',0,'f7',0,'f8',0,'f9',0,'fDynamic',0,'hyper',0);
  for iNode = 2:numNode, namePara(iNode) = namePara(1); numPara(iNode) = numPara(1); end      

  for iNode = 1:numNode-1 
    namePara(iNode).all  = nameList.para;
    numPara(iNode).all   = length(namePara(iNode).all);

    if iNode == 1         
      paraTarget               = paraTargetList(1); 
      namePara(iNode).termName = 'initial';
      numPara(iNode).termName  = 'initial';         
    else
      paraTarget               = paraTargetList(iNode-1);
      namePara(iNode).termName = paraTarget.termName;
      numPara(iNode).termName  = paraTarget.termName;
    end
    
    for name = namePara(iNode).all
      fX = strcat('f',num2str(paraTarget.(name).flag));
      [namePara(iNode),numPara(iNode)] = setForNumName(fX,name,namePara(iNode),numPara(iNode));
      if paraTarget.(name).flag ~= 0, [namePara(iNode),numPara(iNode)] = setForNumName('fDynamic',name,namePara(iNode),numPara(iNode)); end
      if useHyperParameter && paraTarget.(name).setHyper,  [namePara(iNode),numPara(iNode)] = setForNumName('hyper',name,namePara(iNode),numPara(iNode)); end %2022.5.30
    end   
  end
  namePara(end) = namePara(end-1); namePara(end).termName = 'final';
  numPara(end)  = numPara(end-1);  numPara(end).termName  = 'final';

  function [namePara,numPara] = setForNumName(fX,name,namePara,numPara)
    numPara.(fX) = numPara.(fX) + 1;
    namePara.(fX)(numPara.(fX)) = name;
  end
end
% Handles emissions in total CH4, 13CH4, CH3D, and 14CH4.
% Ryo Fujita 2024

classdef Emission 
  properties
    tspan 
    CH4  = SourceWorking;
    C13  = SourceWorking;
    CHD  = SourceWorking;
    C14  = SourceWorking;
    CH4_org
    d13C  
    dD    
    D14C  
  end

  methods
    function obj = SetTspan(obj,inputTspan)
      obj.tspan = inputTspan;
    end

%% 
    function obj = Emission(model,varargin) 
      if nargin == 0, return; end 

      if model.isPriorOnly
        obj.tspan = model.para(1,1).tspan;
        obj = obj.SetEmission(model,model.para(1,1),varargin{:});
        return;
      end

      obj(size(model.para)) = obj; 
      for iRun = 1:size(obj,1)  
        for iNode = 1:size(obj,2)
          model.currentRun   = iRun; 
          model.currentNode  = iNode;  
          obj(iRun,iNode).tspan = model.para(iRun,iNode).tspan; 

          if iNode == 1 || model.isReconstruction
            obj(iRun,iNode) = obj(iRun,iNode).SetEmission(...
              model,model.para(iRun,iNode),varargin{:}); 
          end

        end
      end
    end

%% 
    function obj = SetEmission(obj,model,para,varargin) 
      if isempty(para.fbb)
        disp('para is empty! SetEmission cannot be performed...'); 
        return; 
      end 

      if model.numData ~= size(para.fbb,1)
        model.numData = size(para.fbb,1); 
      end
 
      if model.numCase ~= size(para.fbb,2)
        model.numCase = size(para.fbb,2); 
      end 
     
      if isempty(varargin)
        categoryList = "tot"; 
      else
        categoryList = varargin{1}; 
      end 
      
      obj = obj.SetEmsINodeCH4(model,para,categoryList);
      obj = obj.SetEmsINodeIsotope("d13C",model,para,categoryList);
      obj = obj.SetEmsINodeIsotope("dD",  model,para,categoryList);
      obj = obj.SetEmsINodeRadiocarbon(model,para,categoryList);
    end

%% varargin: categoryList
    function obj = SetEmsINodeCH4(obj,model,para,varargin) 
      model.para = []; %for saving memory
      numCase  = model.numCase;
      numData  = model.numData; 
      paraName = string(fieldnames(para));

      obj.CH4_org   = Emission.ReadFileCH4(model,para.tspan);
      obj.CH4.tspan = para.tspan;

      for cname = string(fieldnames(SourcePrimary)).'
        obj.CH4.(cname) = zeros(numData,numCase);
        switch cname
          case {'agr', 'rumi', 'rice', 'wast'}
            if ismember('fbio',paraName) && ~isempty(para.fbio)
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*para.fbio;
            elseif ismember('fanth_bio',paraName) && ~isempty(para.fanth_bio)
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*para.fanth_bio;
            else
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*ones(1,numCase);
            end
            
          case {'wet', 'trmt', 'anim'}
            if ismember('fbio',paraName) && ~isempty(para.fbio) 
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*para.fbio;
            elseif ismember('fnatr_bio',paraName) && ~isempty(para.fnatr_bio)
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*para.fnatr_bio;
            else
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*ones(1,numCase);
            end
            
          case {'energy', 'gas', 'coal', 'rco', 'otherff'}
            if sum(ismember(paraName,'fff'),1) == 1 && ~isempty(para.fff)
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*para.fff;
            elseif  sum(ismember(paraName,'fanth_ff'),1) == 1 && ~isempty(para.fanth_ff)
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*para.fanth_ff;
            else
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*ones(1,numCase);
            end
            
          case 'bb'
            if ismember('fbb',paraName) && ~isempty(para.fbb)
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*para.fbb;
            else
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*ones(1,numCase);
            end
            
          case 'geo'
            if ismember('Egeo',paraName) && ~isempty(para.Egeo)
              obj.CH4.(cname)  =  ones(numData,1) .* para.Egeo;
            else
              obj.CH4.(cname)  =  obj.CH4_org.(cname)  .*ones(1,numCase);
            end
        end
      end

      obj.CH4 = obj.CH4.computeEmsSecondary(varargin{:});
      obj.CH4 = obj.CH4.computeEmsWorking;   
    end

%% 
    function y = getCH4only(obj,model,varargin)
      if ~isempty(varargin)
        categoryList = varargin{1}; 
      else
        categoryList = []; 
      end

      obj(size(model.para)) = obj;
      for iRun = 1:size(obj,1) 
        for iNode = 1:size(obj,2)
          if isempty(model.para(iRun,iNode).fbb)
            continue; 
          end
 
          if model.numCase ~= size(model.para(iRun,iNode).fbb,2)
            model.numCase = size(model.para(iRun,iNode).fbb,2); 
          end

          obj(iRun,iNode).tspan = model.para(iRun,iNode).tspan;
          obj(iRun,iNode) = obj(iRun,iNode).SetEmsINodeCH4(...
            model,model.para(iRun,iNode),categoryList);
        end
      end

      y = SourceWorking(obj,"CH4");
    end

%%
    function obj = SetEmsINodeIsotope(obj,gasIso,model,para,categoryList)
      switch gasIso
        case "d13C"
          gasEms = "C13"; Rstd = MassConstants.R13C_std;
        case "dD"
          gasEms = "CHD"; Rstd = MassConstants.RD_std;
      end

      for name = categoryList
        cname    = Isotopes.cnameUse(name,model,para);
        namePara = strcat(gasIso,cname);
        switch name
          case 'tot'
            obj.(gasIso).(name) = SourceWorking.(namePara)(para,obj.CH4);

          case 'bio' 
            if ismember(namePara,model.nameList.para)
              obj.(gasIso).(name) = para.(namePara);
            elseif ~isempty(obj.(gasIso))
              fieldnameList = string(fieldnames(obj.(gasIso))).'; 
              if ~ismember("anth_bio",fieldnameList) || isempty(obj.(gasEms).anth_bio)
                obj = obj.SetEmsINodeIsotope(gasIso,model,para,"anth_bio"); 
              end

              if ~ismember("natr_bio",fieldnameList) || isempty(obj.(gasEms).natr_bio)
                obj = obj.SetEmsINodeIsotope(gasIso,model,para,"natr_bio"); 
              end

              obj.(gasEms).bio = obj.(gasEms).anth_bio + obj.(gasEms).natr_bio; 
              obj.(gasIso).bio = (obj.CH4.anth_bio.*obj.(gasIso).anth_bio + ...
                obj.CH4.natr_bio.*obj.(gasIso).natr_bio)./obj.CH4.bio; 
              continue 
            elseif ~isempty(para)
              fieldnameList = string(fieldnames(para)).';
              for paraName = fieldnameList
                if contains(paraName,gasIso)
                  srcName = extractAfter(paraName,gasIso);
                  obj.(gasIso).(srcName) = para.(paraName);
                end
              end
              obj = obj.SetEmsINodeIsotope(gasIso,model,para,name);
            else
              error("para is empty!")
            end

          case 'ff'
            if ismember(namePara,model.nameList.para)
              obj.(gasIso).(name) = para.(namePara);
            elseif ~isempty(obj.(gasIso)) 
              fieldnameList = string(fieldnames(obj.(gasIso))).';
              if ~ismember("anth_ff",fieldnameList) || isempty(obj.(gasEms).anth_ff)
                obj = obj.SetEmsINodeIsotope(gasIso,model,para,"anth_ff"); 
              end

              if ~ismember("geo",fieldnameList) || isempty(obj.(gasEms).geo)
                obj = obj.SetEmsINodeIsotope(gasIso,model,para,"geo"); 
              end

              obj.(gasEms).ff = obj.(gasEms).anth_ff + obj.(gasEms).geo;
              obj.(gasIso).ff = (obj.CH4.anth_ff.*obj.(gasIso).anth_ff + ...
                obj.CH4.geo.*obj.(gasIso).geo)./obj.CH4.ff; 
              continue
            elseif ~isempty(para)
              fieldnameList = string(fieldnames(para)).';
              for paraName = fieldnameList
                if contains(paraName,gasIso)
                  srcName = extractAfter(paraName,gasIso);
                  obj.(gasIso).(srcName) = para.(paraName);
                end
              end
              obj = obj.SetEmsINodeIsotope(gasIso,model,para,name);
            else
              error("para is empty!")
            end

          case 'geo'
            if isempty(obj.CH4.(name)), obj.CH4.(name) = obj.CH4.natr_ff; end
            obj.(gasIso).(name) = para.(namePara);

          case 'natr_ff'
            namePara = strrep(namePara,"natr_ff","geo");
            obj.(gasIso).(name) = para.(namePara);

          case 'npr'
            obj.(gasEms).(name) = zeros(size(obj.CH4.(name))); continue 

          otherwise
            obj.(gasIso).(name) = para.(namePara);
        end

        if isempty(obj.(gasIso).(name))
          error(strcat("obj.",gasIso,".",name," is empty!")); 
        end

        obj.(gasEms).(name) = obj.CH4.(name) .* Rstd .* ...
          (1 + obj.(gasIso).(name)./1000) ./(Rstd.*(1 + obj.(gasIso).(name)./1000) + 1);
      end
    end

%% 
    function obj = SetEmsINodeRadiocarbon(obj,model,para,categoryList) 
      categoryList = convertCatListRadiocarbon(para,categoryList); 
      model.para = []; %clear memory
      numData  = model.numData; 
      numCase  = model.numCase; 

      for cname = categoryList
        obj.C14.(cname)  = zeros(numData ,numCase);
        obj.D14C.(cname) = zeros(numData ,numCase);

        switch cname
          case 'tot'
            obj.C14.tot = obj.getC14tot; 
          case 'npr'
            obj.C14.npr = Isotopes.ComputeEmissionNPR(model,para); 
          case {'anth_ff','geo','ff'} 
            obj.D14C.(cname)(:,:) = -1000;
          case {'anth_bio','natr_bio','bb','bio'}
            fieldnameList = string(fieldnames(obj.d13C)).';
            if ismember(cname,fieldnameList)
              d13Cs = obj.d13C.(cname);
            elseif ~isempty(para.(strcat("d13C",cname)))
              d13Cs = para.(strcat("d13C",cname));
            else
              error(strcat("Empty ",cname));
            end

            [obj.C14.(cname),obj.D14C.(cname)] = Isotopes.ReadIsotopeBioWithTau(...
              "C14",obj.CH4.(cname),d13Cs,model,para,cname); 
          otherwise
            error("Not yet prepared!")
        end
      end  
    end

%%
    function y = getC14tot(obj)
      if ~isempty(obj.C14.anth_bio)
        y = obj.C14.anth_ff + obj. C14.geo +  ...
          obj.C14.anth_bio + obj.C14.natr_bio + obj.C14.bb + obj.C14.npr; 
      elseif ~isempty(obj.C14.bio)
        y = obj.C14.anth_ff + obj. C14.geo + ...
         obj.C14.bio + obj.C14.bb + obj.C14.npr; 
      else
        error("isempty(obj.C14.anth_bio) && isempty(obj.C14.bio)")
      end
    end

%%
    function y = isFinite(obj)
      y = all(isfinite(obj.CH4.tot),1) & all(isfinite(obj.C13.tot),1) & ...
        all(isfinite(obj.CHD.tot),1) & all(isfinite(obj.C14.tot),1);
    end

%%
    function obj = SetEmissionLassey07(obj,model)
      iNode = model.currentNode;
      ECH4_ave = 560; ECH4_std = 40; 
      FracFF_ave = 30; FracFF_std = 2.3;
      
      obj.d13C.bio = -61 .* ones(model.numData,model.numCase);

      if iNode == 1
        obj.CH4.tot  = (ECH4_ave + ECH4_std*randn(1,model.numCase))...
          .* ones(model.numData,1);
        obj.CH4.ff   = obj.CH4.tot .* ((FracFF_ave + ...
          FracFF_std * randn(1,model.numCase))/100.* ones(model.numData,1));
        obj.CH4.bio  = obj.CH4.tot - obj.CH4.ff;

        R_std        = MassConstants.R13C_std; 
        obj.C13.ff   = obj.CH4.ff .* R_std .* ...
          (1 + -45./1000) ./(R_std.*(1 -45./1000) + 1);
        obj.C13.bio  = obj.CH4.bio .* R_std .* ...
          (1 + -61./1000) ./(R_std.*(1 -61./1000) + 1);
        obj.C13.tot  = obj.C13.bio + obj.C13.ff;

        R_std        = MassConstants.RD_std; 
        obj.CHD.ff   = obj.CH4.tot .* R_std .* ...
          (1 + -200./1000) ./(R_std.*(1 -200./1000) + 1); 
        obj.CHD.bio  = obj.CH4.tot .* R_std .* ...
          (1 + -260./1000) ./(R_std.*(1 -260./1000) + 1); 
        obj.C13.tot  = obj.C13.bio + obj.C13.ff;
      else
        judge_oldIsMatchAll =  model.judge(iNode-1).isMatch.all;
        if model.judge(iNode-1).numMatch.all == 0
          judge_oldIsMatchAll = true(1,length(judge_oldIsMatchAll)); 
        end

        iLast = size(model.ems(iNode-1).CH4.tot,1);
        for gasID = ["CH4","C13","CHD"]
          obj.(gasID).bio = model.ems(iNode-1).(gasID).bio(iLast,judge_oldIsMatchAll) .* ones(model.numData,1);
          obj.(gasID).ff  = model.ems(iNode-1).(gasID).ff(iLast,judge_oldIsMatchAll) .* ones(model.numData,1);
          obj.(gasID).tot = obj.(gasID).bio + obj.(gasID).ff;
        end
      end
      obj.C14.tot  = obj.SetEmissionC14Lassey07(model);
    end

%%
    function y = SetEmissionC14Lassey07(obj,model)
        [obj.C14.bio,obj.D14C.bio] = Isotopes.ReadIsotopeBioWithTau(...
          "C14",obj.CH4.bio,obj.d13C.bio,model,"bio"); 
        obj.C14.npr = Isotopes.ComputeEmissionNPR(model); 
        y = obj.C14.bio + obj.C14.npr;
    end
  end

  methods(Static) 
    function CH4_org = ReadFileCH4(model,tspan)
      CH4_org = Emission.readPriorInvFile(model,tspan);
      CH4_org = Emission.setPriorNatr(model,CH4_org);                   
      CH4_org = Emission.setPriorBB(model,CH4_org);                       
      CH4_org = SourceSecondary.computeEmission(CH4_org);
      CH4_org = SourceWorking.computeEmission(CH4_org);    
    end

%% 
    function CH4 = readPriorInvFile(model,varargin)
      [data,header] = MyFunc.readText(model.emsCH4File,10,true);    
      CH4_all.time = data(:,1);

      if ~isempty(varargin)
        CH4.tspan = varargin{1}.'; 
      else
        CH4.tspan = CH4_all.time.'; 
      end

      for i = 2:size(data,2) 
        cname = header(i);
        CH4_all.(cname) = data(:,i); 
        CH4.(cname) = Emission.interp(model,CH4_all.time,CH4_all.(cname),CH4.tspan);
      end 
    end

%% 
    function y = interp(model,tspan_org,ECH4_org,tspan) 
      t1 = tspan <= model.yrInitial;
      t2 = model.yrInitial < tspan & tspan <= tspan_org(end);
      t3 = tspan > tspan_org(end);
      i1 = find(t1);
      i2 = find(t2);
      i3 = find(t3);

      if sum(i1) >= 1
        y(i1,:) = interp1(tspan_org, ECH4_org, model.yrInitial, 'next') ...
          .*ones(size(i1,2),1); 
      end 

      if sum(i2) >= 1 
        if isempty(i1)
          y = interp1(tspan_org, ECH4_org, tspan(i2), 'linear');
        else
          y(i2,:) = interp1(tspan_org, ECH4_org, tspan(i2), 'linear');
        end
      end

      if sum(i3) >= 1
        y(i3,:) = ECH4_org(end).*ones(size(i3,2),1); 
      end

      if contains(model.Info.inventory,"constAft")
        yrBase = str2double(extractAfter(model.Info.inventory,"constAft"));
        logic = tspan > yrBase+0.5;
        if any(logic)
          y(i3,:) = interp1(tspan_org, ECH4_org, yrBase+0.5, 'linear'); 
        end
      end
    end

%%
    function CH4 = setPriorNatr(model,CH4)
      numData = length(CH4.tspan);

      if size(CH4.tspan,1) == 1
        tspan_wet = CH4.tspan.'; 
      else
        tspan_wet = CH4.tspan; 
      end

      switch model.emsNatr
        case {'S20',''} 
          CH4.wet(1:numData,1) = 306;
          CH4.anim(1:numData,1) = 2;
          CH4.trmt(1:numData,1)= 9;
          CH4.geo(1:numData,1) = 40;
        case {'S20_wet147to180'}
          % 159 Tg/yr Freshwater, 147 Tg/yr wetland in Saunois et al. 2020, Table 3 
          % for 2000-2009
          CH4.wet = 159 + ...
            interp1([1650,1850,2015.5],[147,147,180],tspan_wet,'linear','extrap');
          CH4.anim(1:numData,1) = 2;
          CH4.trmt(1:numData,1)= 9;
          CH4.geo(1:numData,1)= 40;
        case {'S20_wetVISIT-Cao'} 
          visit = readtable("emission/VISIT/visit_ver.2021.1_CH4Wetl_Cao_annual.csv");
          tspan = [1650;1900.5;visit.tspan];
          ems = [166.0;166.0;visit.ems]; %average for 1901-1910
          CH4.wet = 159 + interp1(tspan,ems,tspan_wet,'linear','extrap');
          CH4.anim(1:numData,1) = 2;
          CH4.trmt(1:numData,1) = 9;
          CH4.geo(1:numData,1) = 40;
        case {'S20_wetVISIT-WH'}
          visit = readtable("emission/VISIT/visit_ver.2021.1_CH4Wetl_WH_annual.csv");
          tspan = [1650;1900.5;visit.tspan];
          ems = [132.0;132.0;visit.ems]; %average for 1901-1910
          CH4.wet = 159 + interp1(tspan,ems,tspan_wet,'linear','extrap');
          CH4.anim(1:numData,1) = 2;
          CH4.trmt(1:numData,1) = 9;
          CH4.geo(1:numData,1) = 40;
        case {'S20_noFreshwater'}
          CH4.wet(1:numData,1)  = 306-159;
          CH4.anim(1:numData,1) = 2;
          CH4.trmt(1:numData,1) = 9;
          CH4.geo(1:numData,1)  = 40;
        case 'L07'
          CH4.wet(1:numData,1)  = 163;
          CH4.anim(1:numData,1) = 20;
          CH4.trmt(1:numData,1) = 15;
          CH4.geo(1:numData,1)  = 19;
        case 'dummy'
          CH4.wet(1:numData,1)  = 0;
          CH4.anim(1:numData,1) = 0;
          CH4.trmt(1:numData,1) = 0;
          CH4.geo(1:numData,1)  = 0;
        otherwise
          error(strcat("No such emsNatr exists!: ",model.emsNatr)) 
      end
    end

%%
    function CH4 = setPriorBB(model,CH4)
      numData = length(CH4.tspan);
      switch model.emsBB
        case {'M17','','M17plusRCO'} %van Marle et al. 2017
          [year_bb, Ebb] = MyFunc.readTextFreefmt("CH4_em_biomass_input4MIPs_1700-2015.txt", '%n %n','\t');
          CH4.bb = Emission.interp(model, year_bb, Ebb, CH4.tspan); 

          if strcmp(model.emsBB,'M17plusRCO')
            CH4.bb = CH4.bb + CH4.rco;
            CH4.rco = zeros(size(CH4.rco));
          end
        case 'S16' %Schwietzke et al., 2016
          CH4.bb(1:numData,1) = 43;
        case 'S20' %Saunois et al., 2020, bottom-up mean for 2000-2009
          CH4.bb(1:numData,1) = 31; 
        case 'dummy'
          CH4.bb(1:numData,1) = 30; 
      end
    end

%% 
    function y = shortName()
      y = "ems";
    end
  end

end

%%
function y = convertCatListRadiocarbon(para,categoryList)
  if ~isempty(para.d13Canth_bio)
    nameList = ["anth_bio","natr_bio","anth_ff","geo","bb","npr"];   
  elseif ~isempty(para.d13Cbio)
    nameList = ["bio","anth_ff","geo","bb","npr"];
  else
    error("isempty(obj.C14.anth_bio) && isempty(obj.C14.bio)")
  end 

  pos = find(ismember(categoryList,"tot")); 
  if isempty(pos)
    y = categoryList; 
    return
  elseif pos ~= length(categoryList) 
    y = [categoryList(~ismember(categoryList,"tot")),categoryList(pos)]; %"tot" should be placed at the last
  else
    y = categoryList;
  end

  if sum(ismember(y,nameList)) ~= length(nameList)
    for name = nameList
      if ~ismember(name,y)
        y = [name,y];
      end
    end
  end
end
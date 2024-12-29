classdef SourceWorking < SourceSecondary
  
  properties
    tspan 
    bio
    ff
    %bb %already defined in SourcePrimary
    npr
    tot
  end

  methods
    function obj = SourceWorking(varargin) 
      if nargin ~= 0 
        emsObj  = varargin{1};
        gasID   = varargin{2};

        if nargin == 3 && ~isempty(varargin{3})
          structList = varargin{3};
        else
          structList = string(fieldnames(emsObj(1,1).(gasID))).';
        end
        clear varargin

        numRun  = size(emsObj,1);
        numNode = size(emsObj,2);
        obj(numRun,numNode) = obj;

        for iRun = 1:numRun
        for iNode = 1:numNode
          obj(iRun,iNode).tspan = emsObj(iRun,iNode).tspan;
          if isempty(emsObj(iRun,iNode).(gasID)), continue; end 
          for name = structList
            if sum(ismember(string(fieldnames(emsObj(iRun,iNode).(gasID))).',name)) == 0
              continue; 
            end
            obj(iRun,iNode).(name) = emsObj(iRun,iNode).(gasID).(name); 
            emsObj(iRun,iNode).(gasID).(name) = [];
          end
        end
        end
      end
    end

%%
    function obj = computeEmsSecondary(obj,varargin)
      obj.anth_ff  = obj.gas + obj.coal + obj.rco  + obj.otherff;
      obj.anth_bio = obj.rumi + obj.rice + obj.wast ;
      obj.natr_ff  = obj.geo; 
      obj.natr_bio = obj.wet  + obj.trmt + obj.anim; 
      obj.anth     = obj.anth_ff + obj.anth_bio + obj.bb;
      obj.natr     = obj.natr_ff + obj.natr_bio;

      if isempty(varargin)
        categoryList = []; 
      else
        categoryList = varargin{1}; 
      end
 
      obj = obj.clearMemory(SourcePrimary,categoryList,["bb","geo"]);
    end

%%
    function obj = computeEmsWorking(obj)
      obj.ff   =  obj.anth_ff + obj.natr_ff;
      obj.bio  =  obj.anth_bio + obj.natr_bio;
      obj.npr  =  zeros(size(obj.bio)); 
      obj.tot  =  obj.anth + obj.natr; 
    end

%% 
    function obj = clearMemory(obj,clearClass,nameList,varargin)
      if ~isempty(varargin)
        saveNameList = varargin{1}; 
      else
        saveNameList = ""; 
      end

      clearNameList = ClearMemory.getFieldList(clearClass,nameList);
      for name = clearNameList
        if ismember(name,saveNameList)
          continue; 
        end

        obj.(name) = []; 
      end 
    end
  end

%%
  methods (Static)
    function inout = computeEmission(inout) 
      inout.ff   =  inout.anth_ff + inout.natr_ff;
      inout.bio  =  inout.anth_bio + inout.natr_bio;
      inout.npr  =  zeros(size(inout.bio));
      inout.tot  =  inout.anth + inout.natr;
    end

%%
    function inout = computeEmissionRadiocarbon(inout) 
      if isempty(inout.natr_ff)
        inout.natr_ff = inout.geo; 
      end

      inout.ff   =  inout.anth_ff + inout.natr_ff;
      inout.bio  =  inout.anth_bio + inout.natr_bio;

      if ~isempty(inout.anth) && ~isempty(inout.natr) && ~isempty(inout.npr)
        inout.tot  =  inout.anth + inout.natr + inout.npr;
      elseif ~isempty(inout.ff) && ~isempty(inout.bb) && ~isempty(inout.bio) ...
      && ~isempty(inout.npr)      
        inout.tot  =  inout.ff + inout.bb + inout.bio + inout.npr;
      else
        error("Please check inout!")
      end
    end

%%  
    function d14Cs = computeDeltaRadiocarbon(E_C14,E_CH4,d13Cs,d14Cs)
      M    = MassConstants;
      Aabs = M.Aabs;
      
      d14Cs.ff  = ((E_C14.ff.* (0.975./(1+d13Cs.ff./1000)).^2 )./E_CH4.ff./Aabs -1).*1000 ;
      d14Cs.bio = ((E_C14.bio.* (0.975./(1+d13Cs.bio./1000)).^2 )./E_CH4.bio./Aabs -1).*1000 ; 
      d14Cs.tot = ((E_C14.bio.* (0.975./(1+d13Cs.bio./1000)).^2 ...
                  + E_C14.bb.* (0.975./(1+d13Cs.bb./1000)).^2 ...
                  + E_C14.npr)./E_CH4.tot./Aabs -1).*1000 ; 
    end 
 
%% 
    function delta = computeStableIsotopeRatio(gasID,E_ISO,E_CH4,delta,varargin)
      switch gasID
        case 'd13C'
          R_std = MassConstants.R13C_std; 
        case 'dD'
          R_std = MassConstants.RD_std; 
      end
      delta.ff  = (E_ISO.ff ./(E_CH4.ff  - E_ISO.ff) ./R_std -1).*1000;
      delta.bb  = (E_ISO.bb ./(E_CH4.bb  - E_ISO.bb) ./R_std -1).*1000;
      delta.bio = (E_ISO.bio./(E_CH4.bio - E_ISO.bio)./R_std -1).*1000;
      delta.tot = (E_ISO.tot./(E_CH4.tot - E_ISO.tot)./R_std -1).*1000;
    end

%% 
    function y = d13Ctot(para,E_CH4)
      if isempty(para.d13Cbio)
        para.d13Cbio = SourceWorking.deltaBIO_wtMean( ...
          E_CH4.anth_bio, para.d13Canth_bio, E_CH4.natr_bio, para.d13Cnatr_bio);
      end

      if isempty(para.d13Cff)
        para.d13Cff  = SourceWorking.deltaFF_wtMean( ...
          E_CH4.anth_ff, para.d13Canth_ff, E_CH4.natr_ff, para.d13Cgeo);
      end

      y  = SourceWorking.deltaTOT_wtMean( ...
        E_CH4.anth_ff+E_CH4.natr_ff, para.d13Cff, E_CH4.anth_bio+E_CH4.natr_bio, ...
        para.d13Cbio, E_CH4.bb, para.d13Cbb);
    end  

%%
    function y = dDtot(para,E_CH4)
      if isempty(para.dDbio)
        para.dDbio = SourceWorking.deltaBIO_wtMean( ...
          E_CH4.anth_bio, para.dDanth_bio, E_CH4.natr_bio, para.dDnatr_bio);
      end

      if isempty(para.dDff)
        para.dDff  = SourceWorking.deltaFF_wtMean( ...
          E_CH4.anth_ff, para.dDanth_ff, E_CH4.natr_ff, para.dDgeo);
      end

      y  = SourceWorking.deltaTOT_wtMean( ...
        E_CH4.anth_ff+E_CH4.natr_ff, para.dDff, E_CH4.anth_bio+E_CH4.natr_bio, ...
        para.dDbio, E_CH4.bb, para.dDbb);
    end 

%%
    function delta = WeightedMeanStableIsotopesDef(objName,E_CH4,delta) 
      delta.ff.def   = SourceWorking.deltaFF_wtMean( ...
        E_CH4.anth_ff, delta.anth_ff.def, E_CH4.natr_ff, delta.natr_ff.def);

      delta.bio.def  = SourceWorking.deltaBIO_wtMean( ...
        E_CH4.anth_bio, delta.anth_bio.def, E_CH4.natr_bio, delta.natr_bio.def);

      delta.tot.def  = SourceWorking.deltaTOT_wtMean( ...
        E_CH4.ff, delta.ff.def, E_CH4.bio, delta.bio.def, E_CH4.bb, delta.bb.def);

      cnameList = ["bio","ff","tot"];
      delta = Isotopes.SetMinMaxDef(objName,delta,cnameList);
    end

%%
    function y = deltaFF_wtMean(E_anthFF,delta_anthFF,E_natrFF,delta_natrFF)
      y = (E_anthFF .* delta_anthFF + E_natrFF .* delta_natrFF)./(E_anthFF + E_natrFF);
    end

%%
    function y = deltaBIO_wtMean(E_anthBIO,delta_anthBIO,E_natrBIO,delta_natrBIO)
      y = (E_anthBIO .* delta_anthBIO + E_natrBIO .* delta_natrBIO)./(E_anthBIO + E_natrBIO);
    end

%%
    function y = deltaTOT_wtMean(E_FF,delta_FF,E_BIO,delta_BIO,E_BB,delta_BB)
      y = (E_FF .* delta_FF + E_BIO .* delta_BIO + E_BB .* delta_BB)./(E_FF + E_BIO + E_BB);
    end

  end
end

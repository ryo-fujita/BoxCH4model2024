classdef SourceSecondary < SourcePrimary
  properties
    anth_ff
    natr_ff 
    anth_bio
    natr_bio
    natr
    anth
  end

  methods (Static)
    function inout = computeEmission(inout) 
      inout.anth_ff        = inout.gas + inout.coal + inout.rco  + inout.otherff;
      
      inout.anth_bio       = inout.rumi + inout.rice + inout.wast ;
      inout.natr_ff        = inout.geo; 
      inout.natr_bio       = inout.wet    + inout.trmt + inout.anim; 
      
      inout.anth           = inout.anth_ff + inout.anth_bio + inout.bb;
      inout.natr           = inout.natr_ff + inout.natr_bio;
    end

%% 
    function inout = clearMemory(inout)
      for name = string(fieldnames(SourceSecondary)).'
        inout.(name) = [];
      end
    end

%%
    function delta = computeStableIsotopeRatio(gasID,E_ISO,E_CH4,delta,varargin)
      switch gasID
        case 'd13C'
          R_std    = MassConstants.R13C_std; 
        case 'dD'
          R_std    = MassConstants.RD_std; 
      end

      delta.anth           = (E_ISO.anth ./(E_CH4.anth  - E_ISO.anth) ./R_std -1).*1000;
      delta.anth_ff        = (E_ISO.anth_ff ./(E_CH4.anth_ff  - E_ISO.anth_ff) ./R_std -1).*1000;
      delta.anth_bio       = (E_ISO.anth_bio./(E_CH4.anth_bio - E_ISO.anth_bio)./R_std -1).*1000;
      delta.natr           = (E_ISO.natr ./(E_CH4.natr  - E_ISO.natr) ./R_std -1).*1000;
      delta.natr_ff        = (E_ISO.natr_ff ./(E_CH4.natr_ff  - E_ISO.natr_ff) ./R_std -1).*1000;
      delta.natr_bio       = (E_ISO.natr_bio./(E_CH4.natr_bio - E_ISO.natr_bio)./R_std -1).*1000;
    end

%%
    function delta = WeightedMeanStableIsotopesDef(objName,E_CH4,delta)  
      delta.anth_ff.def   = (E_CH4.energy .* delta.energy.def + E_CH4.rco .* delta.rco.def + E_CH4.otherff .* delta.otherff.def)./E_CH4.anth_ff;
      delta.anth_bio.def  = (E_CH4.agr .* delta.agr.def + E_CH4.wast .* delta.wast.def)./E_CH4.anth_bio ;
      delta.natr_ff.def   = delta.geo.def;
      delta.natr_bio.def  = (E_CH4.wet .* delta.wet.def + E_CH4.trmt .* delta.trmt.def + E_CH4.anim .* delta.anim.def)./E_CH4.natr_bio;  

      delta.anth.def      = (E_CH4.anth_ff .* delta.anth_ff.def + E_CH4.anth_bio .* delta.anth_bio.def + E_CH4.bb .* delta.bb.def)./E_CH4.anth;
      delta.natr.def      = (E_CH4.natr_ff .* delta.natr_ff.def + E_CH4.natr_bio .* delta.natr_bio.def)./E_CH4.natr;

      cnameList = string(fieldnames(SourceSecondary)).';
      delta = Isotopes.SetMinMaxDef(objName,delta,cnameList);
    end

%% 
    function d14Cs = computeDeltaRadiocarbon(E_C14,E_CH4,d13Cs,d14Cs)
      M    = MassConstants;
      Aabs = M.Aabs;

      d14Cs.anth           = ((E_C14.anth .* (0.975./(1+d13Cs.anth ./1000)).^2 )          ./E_CH4.anth ./Aabs -1).*1000 ;  
      d14Cs.anth_ff        = ((E_C14.anth_ff .* (0.975./(1+d13Cs.anth_ff ./1000)).^2 )    ./E_CH4.anth_ff ./Aabs -1).*1000 ;   
      d14Cs.anth_bio       = ((E_C14.anth_bio.* (0.975./(1+d13Cs.anth_bio./1000)).^2 )    ./E_CH4.anth_bio./Aabs -1).*1000 ;  
      d14Cs.natr           = ((E_C14.natr .* (0.975./(1+d13Cs.natr ./1000)).^2 )          ./E_CH4.natr ./Aabs -1).*1000 ; 
      d14Cs.natr_ff        = ((E_C14.natr_ff .* (0.975./(1+d13Cs.natr_ff./1000)).^2 )     ./E_CH4.natr_ff./Aabs -1).*1000 ;   
      d14Cs.natr_bio       = ((E_C14.natr_bio.* (0.975./(1+d13Cs.natr_bio./1000)).^2 )    ./E_CH4.natr_bio./Aabs -1).*1000 ;   
    end
  end
end

% A class that aggregates the properties of name lists. 
% It is located in constants/.
% Ryo Fujita 2024

classdef MyNameList
  properties (Constant)
    paraBIO_default     = ["fbio","fanth_ff","Egeo","fbb","tau","phi","life","KIEC","KIED"];
    paraBIO_separated   = ["fanth_bio","fnatr_bio","fanth_ff","Egeo","fbb","tau","phi","life","KIEC","KIED"];
    paraBIO_separateWet = ["fanth_bio","fwet","fanth_ff","Egeo","fbb","tau","phi","life","KIEC","KIED"];
    paraBIO_dummy       = ["fanth_bio","fanth_ff","fbb","tau","phi","life","KIEC","KIED"];

    isoBIO_default      = ["d13Cbio","d13Canth_ff","d13Cgeo","d13Cbb","dDbio","dDanth_ff","dDgeo","dDbb"];
    isoBIO_3source      = ["d13Cbio","d13Cff","d13Cbb","dDbio","dDff","dDbb"];
    isoBIO_separated    = ["d13Canth_bio","d13Cnatr_bio","d13Canth_ff","d13Cgeo","d13Cbb","dDanth_bio","dDnatr_bio","dDanth_ff","dDgeo","dDbb"];
    isoBIO_separateWet  = ["d13Canth_bio","d13Cwet","d13Canth_ff","d13Cgeo","d13Cbb","dDanth_bio","dDwet","dDanth_ff","dDgeo","dDbb"];

    emsBIO_default      = ["tot","bio","ff","bb"];
    emsBIO_separated    = ["tot","anth_bio","natr_bio","ff","bb"];

    emsBIO2_default     = ["bio","anth_ff","geo","bb"]; 
    emsBIO2_separated   = ["anth_bio","natr_bio","anth_ff","geo","bb"];

    frac_default        = ["bio","ff","bb"];
    fracBIO_default     = ["bio","ff","bb"];
    fracBIO_separated   = ["anth_bio","natr_bio","anth_ff","geo","bb"];

    %ssp = ["IMAGE-ssp119","IMAGE-ssp126","MESSAGE-GLOBIOM-ssp245","AIM-ssp370","REMIND-MAGPIE-ssp534-over","REMIND-MAGPIE-ssp585"];
  end

  methods(Static)
    function y = para(model)
      nameList = model.nameParaList(model.numNode).fDynamic.';

      if sum(ismember(nameList,'fbio'),'all') == 1 
        y = MyNameList.paraBIO_default;
      elseif sum(ismember(nameList,'fanth_bio'),'all') == 1 && sum(ismember(nameList,'fnatr_bio'),'all') == 1
        y = MyNameList.paraBIO_separated;
      elseif sum(ismember(nameList,'fanth_bio'),'all') == 1 && sum(ismember(nameList,'fwet'),'all') == 1
        y = MyNameList.paraBIO_separateWet;
      elseif  sum(ismember(nameList,'fanth_bio'),'all') == 1 
        y = MyNameList.paraBIO_dummy;
      else
        y = []; disp('Unproper parameter name setting for bio!')
      end

      if ismember("floss",nameList) && ~ismember("life",nameList) %2023.5.16
        y = replace(y,"life","floss");
      end
    end

%%
    function y = iso(model)
      nameList = model.nameParaList(model.numNode).all.';

      if sum(ismember(nameList,'d13Cbio'),'all') == 1 && sum(ismember(nameList,'d13Canth_ff'),'all') == 1
        y = MyNameList.isoBIO_default;
      elseif sum(ismember(nameList,'d13Cbio'),'all') == 1 && sum(ismember(nameList,'d13Cff'),'all') == 1
        y = MyNameList.isoBIO_3source;
      elseif sum(ismember(nameList,'d13Canth_bio'),'all') == 1 && sum(ismember(nameList,'d13Cnatr_bio'),'all') == 1
        y = MyNameList.isoBIO_separated;
      elseif sum(ismember(nameList,'d13Canth_bio'),'all') == 1 && sum(ismember(nameList,'d13Cwet'),'all') == 1
        y = MyNameList.isoBIO_separateWet;
      else
        y = []; disp('Unproper parameter name setting for bio!')
      end
    end

%%
    function y = ems(model)
      nameList = model.nameParaList(model.numNode).fDynamic.';

      if sum(ismember(nameList,'fbio'),'all') == 1 
        y = MyNameList.emsBIO_default;
      elseif sum(ismember(nameList,'fanth_bio'),'all') == 1 && sum(ismember(nameList,'fnatr_bio'),'all') == 1
        y = MyNameList.emsBIO_separated;
      elseif sum(ismember(nameList,'fanth_bio'),'all') == 1 && sum(ismember(nameList,'fwet'),'all') == 1
        y = MyNameList.emsBIO_separateWet;
      elseif  sum(ismember(nameList,'fanth_bio'),'all') == 1 
        y = MyNameList.emsBIO_default;
      else
        y = []; disp('Unproper parameter name setting for bio!')
      end
    end

%%
    function y = emsSec(model)
      nameList = model.nameParaList(model.numNode).fDynamic.';

      if sum(ismember(nameList,'fbio'),'all') == 1 
        y = MyNameList.emsBIO2_default;
      elseif sum(ismember(nameList,'fanth_bio'),'all') == 1 && sum(ismember(nameList,'fnatr_bio'),'all') == 1
        y = MyNameList.emsBIO2_separated;
      else
        y = []; disp('Unproper parameter name setting for bio!')
      end
    end

%%
    function y = frac(model)
      nameList = model.nameParaList(model.numNode).fDynamic.';

      if sum(ismember(nameList,'fbio'),'all') == 1 
        y = MyNameList.fracBIO_default;
      elseif sum(ismember(nameList,'fanth_bio'),'all') == 1 && sum(ismember(nameList,'fnatr_bio'),'all') == 1
        y = MyNameList.fracBIO_separated;
      elseif sum(ismember(nameList,'fanth_bio'),'all') == 1 && sum(ismember(nameList,'fwet'),'all') == 1
        y = MyNameList.fracBIO_separateWet;
      elseif  sum(ismember(nameList,'fanth_bio'),'all') == 1 
        y = MyNameList.fracBIO_default;
      else
        error('Unproper parameter name setting for bio!')
      end
    end
  end
end


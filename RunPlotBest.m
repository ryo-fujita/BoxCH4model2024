%
classdef RunPlotBest
  properties 
    workdirList
    lgdNameList
    model
  end
  
  methods
    function obj = RunPlotBest(varargin)
      addpath constants
      
      if ~isempty(varargin)
        obj.workdirList = varargin{1};
      else
        obj.workdirList = ".";   
      end
      obj.model = Model(obj.workdirList(1),"Posterior",9);
      obj.lgdNameList = Plot.setlgdNameList(obj.workdirList,obj.model);
    end

%%
    function Prctile(obj)
      plt = Prctile(obj.workdirList,obj.model,"Atmosphere",["CH4","d13C","dD","D14C"]);
      plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1740 2020]}); close(gcf); 
      plt.plot('periodList',{[1980 2016]}); close(gcf);
      clear plt

      plt = Prctile(obj.workdirList,obj.model,"Parameter",obj.model.nameList.para);
      if ~isempty(MyNameList.para(obj.model))
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1850 2020]}, ...
          'categoryList',MyNameList.para(obj.model)); close(gcf); 
      end
      if ~isempty(MyNameList.iso(obj.model))
        plt.plot('categoryList',MyNameList.iso(obj.model)); close(gcf); 
      end
      if ismember("floss",obj.model.nameList.para)
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1850 2020]}, ...
          'categoryList',["floss","life"]); close(gcf); 
      end
      clear plt

      plt = Prctile(obj.workdirList,obj.model,"ECH4", ...
       ["tot","bio","ff","bb","anth_ff","geo","bb","anth_bio","natr_bio",...
        "rice","rumi","wast","coal","gas"]);
      plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1850 2020]}, ...
        'categoryList',["tot","bio","anth_ff","geo","bb"]); close(gcf);
      plt.plot('categoryList',["tot","bio","ff","bb"]); close(gcf);
      plt.plot('categoryList',["tot","anth_bio","anth_ff"]); close(gcf);
      plt.plot('categoryList',["tot","natr_bio","geo"]); close(gcf);
      plt.plot('categoryList',["tot","anth_bio","natr_bio","anth_ff","bb"]); close(gcf);
      plt.plot('categoryList',["rice","rumi","wast","coal","gas","bb"],'subPlotLayout',[2,3]); 
      close(gcf); clear plt
      
      plt = Prctile(obj.workdirList,obj.model,"SourceFraction",["bio","ff","bb"]);
      plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1850 2020]}); 
      close(gcf); clear plt
            
      plt = Prctile(obj.workdirList,obj.model,"d13CdDKIECKIEDtot",["d13Ctot","KIEC","dDtot","KIED"]);
      plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1850 2020]}); 
      close(gcf); clear plt
    end

%% 
    function ParaAllSuppliment(obj)
      categoryList = obj.model.nameList.para;

      %% create Fig.2
      plt = Prctile(obj.workdirList,obj.model,"Parameter",categoryList);
      fig = plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1850 2020]}, ...
       'statsNameList',["prc68","ave"],'subPlotLayout',[4 5],'setSaveFig','no');
      ax = findobj(gcf,'type','ax');

      for i = 1:length(ax)
        plotCharacter(length(ax)-i+1,ax(i),10,'Prctile');
        ax(i).FontSize = 9; ax(i).YLabel.FontSize = 10; ax(i).XLabel.FontSize = 10;
        if i ~= 3, ax(i).XLabel.String = ''; end
        if i == 3, ax(i).XLabel.FontSize = 13; end
      end

      fig.Units = 'centimeter'; 
      fig.Position(3) = fig.PaperSize(1); 
      fig.Position(4) = fig.PaperSize(2); 
      pause(0.5);

      for i = 1:length(ax), ax(i).Position(4) = 0.12; end
      lgd = findobj(gcf,'type','legend');  
      lgd.Position(1) = 0.4;

      saveas(fig,'Fig2_paraAllTimeseries.png');
      savefig(fig,'Fig2_paraAllTimeseries.fig');

      %% create Fig.3
      plt = Histogram(obj.workdirList,obj.model,"Parameter",categoryList);
      fig = plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[2003 2012]}, ...
        'subPlotLayout',[4 5],'setSaveFig','no');
      ax = findobj(gcf,'type','ax');
      fontSize = 11;

      fig.Units = 'centimeter'; 
      fig.Position(3) = fig.PaperSize(1); 
      fig.Position(4) = fig.PaperSize(2); 
      pause(0.5);

      for i = 1:length(ax)
        plotCharacter(length(ax)-i+1,ax(i),fontSize,'Histogram');
        ax(i).FontSize = fontSize-1; 
        ax(i).YLabel.FontSize = fontSize;
        ax(i).Title.FontSize = fontSize ;

        ax(i).Position(3) = 0.11;
        ax(i).Position(4) = 0.14; 
        pause(0.2);

        ax(i).Title.Position(2) = ax(i).YLim(2) + ax(i).YLim(2)/40; 

        if i > 5, ax(i).XLabel.String = ''; end

        if i/5 <= 1
          ax(i).Position(2) = ax(i).Position(2) - 0.03;
        elseif i/5 <= 2
          ax(i).Position(2) = ax(i).Position(2) - 0.02;
        elseif i/5 <= 3
          ax(i).Position(2) = ax(i).Position(2) - 0.01;
        end
      end
      an = annotation('textarrow',[0.08 0.08],[0.6 0.6],'String','Probability', ...
        'FontSize',14,'TextRotation',90,'LineStyle','none','HeadStyle','none');

      obj.addPriorUniform(obj.model,categoryList);
      saveas(fig,'Fig3_paraAllHistogram.png');
      savefig(fig,'Fig3_paraAllHistogram.fig');
    end

%% 
    function Histogram(obj)
      periodListCell = [{[1700 1750; 1750 1900; 1900 1950; 1950 1960; 1960 1970; ...
        1970 1980; 1980 1990; 1990 1995; 1995 2000; 2000 2005; 2005 2010; 2010 2015]}, ...
        {[1940 1950; 1986 2000; 2003 2012]}];
      
      for periodList = periodListCell
        plt = Histogram(obj.workdirList,obj.model,"Parameter",MyNameList.para(obj.model));
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); 
        close(gcf); clear plt
        
        plt = Histogram(obj.workdirList,obj.model,"Parameter",MyNameList.iso(obj.model));
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); 
        close(gcf); clear plt
        
        plt = Histogram(obj.workdirList,obj.model,"ECH4",["tot","anth_bio","natr_bio","anth_ff","geo","bb"]);
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); 
        close(gcf); clear plt
        
        plt = Histogram(obj.workdirList,obj.model,"SourceFraction",["bio","ff","bb"]);
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); 
        close(gcf); clear plt
      end
    end

%%
    function ParaTimeVar(obj)
      if ~obj.model.useParaTimeVar, return; end

      obj.model.lagFilterTerm = obj.model.numNode;
      para = Smoothing.loadPosterior("para",obj.model);
      plt  = Prctile(obj.workdirList,obj.model,"ParaTimeVar",obj.model.nameList.para, ...
        'plotObj',para.IAVpercent); 
      if ~isempty(MyNameList.para(obj.model))
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1740 2020]}, ...
        'categoryList',MyNameList.para(obj.model)); 
        close(gcf); 
      end

      plt = Prctile(obj.workdirList,obj.model,"ParaTimeVar",obj.model.nameList.para,...
       'plotObj',para.IAVpercent); 
      if ~isempty(MyNameList.iso(obj.model))
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1740 2020]},...
        'categoryList',MyNameList.iso(obj.model)); 
        close(gcf); 
      end
    end

%% 
    function ParaHyper(obj)
      if ~obj.model.useHyperParameter, return; end
      obj.model.lagFilterTerm = obj.model.numNode;
      para = Smoothing.loadPosterior("para",obj.model);

      plt = Prctile(obj.workdirList,obj.model,"HyperParameter",obj.model.nameList.para,...
       'plotObj',para.hyper); 

      if ~isempty(MyNameList.para(obj.model))
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1740 2020]},...
        'categoryList',MyNameList.para(obj.model)); 
        close(gcf); 
      end

      if ~isempty(MyNameList.iso(obj.model))
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',{[1740 2020]},...
        'categoryList',MyNameList.iso(obj.model)); 
        close(gcf); 
      end
    end

%% 
    function BarPlot(obj)
      plt = BarPlot(obj.workdirList,obj.model, ...
        "ECH4Life",["tot","anth_bio","natr_bio","anth_ff","bb","life"]);
      plt.plot('plotPriorStyle','yes','lgdNameList',obj.lgdNameList,...
        'calcType',"PeriodDiff",'periodList',{[1991 1998; 1999 2006]}); 
      plt.plot('periodList',{[1999 2006; 2007 2014]}); clear plt
        
      plt = BarPlot(obj.workdirList,obj.model,"ECH4Life",["tot","bio","ff","bb","life"]);
      plt.plot('plotPriorStyle','yes','lgdNameList',obj.lgdNameList, ...
        'calcType',"PeriodDiff",'periodList',{[1991 1998; 1999 2006]}); 
      plt.plot('periodList',{[1999 2006; 2007 2014]}); clear plt

      if length(obj.workdirList) == 1
        periodListCell = {[1986 2000; 2003 2012]};
      else
        periodListCell = [{[1986 2000]},{[2003 2012]}];
      end
    
      for periodList = periodListCell
        plt = BarPlot(obj.workdirList,obj.model,"Parameter",MyNameList.para(obj.model));
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); clear plt
        
        plt = BarPlot(obj.workdirList,obj.model,"Parameter",MyNameList.iso(obj.model));
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); clear plt
        
        plt = BarPlot(obj.workdirList,obj.model,...
         "ECH4",["tot","anth_bio","natr_bio","anth_ff","geo","bb"]);
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); clear plt
        
        plt = BarPlot(obj.workdirList,obj.model,"SourceFraction",["bio","ff","bb"]);
        plt.plot('lgdNameList',obj.lgdNameList,'periodList',periodList); clear plt
      end
    end

%%
    function Gplotmatrix(obj)
      if length(obj.workdirList) == 1
        plotNm = "plot"; 
      else
        plotNm = "plotMean"; 
      end

      plt = Gplotmatrix(obj.workdirList,obj.model,"ECH4",["tot","anth_bio","natr_bio","anth_ff"]);
      plt.(plotNm)('lgdNameList',obj.lgdNameList,'periodList',{[2003 2012]},...
        'categoryList',["tot","anth_bio","anth_ff"]); close(gcf);
      plt.(plotNm)('categoryList',["tot","anth_bio","natr_bio"]); 
      close(gcf);

      plt.(plotNm)('periodList',{[1986 2000]},...
        'categoryList',["tot","anth_bio","anth_ff"]); 
      close(gcf);

      plt.(plotNm)('categoryList',["tot","anth_bio","natr_bio"]); 
      close(gcf); 
      clear plt

      plt = Gplotmatrix(obj.workdirList,obj.model,"Atmosphere",["CH4","d13C","dD","D14C"]);
      plt.(plotNm)('lgdNameList',obj.lgdNameList,'calcType',"LinearTrendDiff",...
       'periodList',{[1991 1998; 1999 2006]});
      close(gcf)

      plt.(plotNm)('periodList',{[1999 2006; 2007 2014]});
      close(gcf); 
      clear plt

      plt = Gplotmatrix(obj.workdirList,obj.model,...
        "ECH4Life",["tot","anth_bio","anth_ff","bb","life"]);
      plt.cutUnit = []; plt.(plotNm)('lgdNameList',obj.lgdNameList,...
        'calcType',"PeriodDiff",'periodList',{[1991 1998; 1999 2006]}); 
      close(gcf);
      plt.(plotNm)('periodList',{[1999 2006; 2007 2014]}); 
      close(gcf); 
      clear plt

      plt = Gplotmatrix(obj.workdirList,obj.model,...
        "Parameter",["fbb","fanth_bio","fanth_ff","tau","phi","life","KIEC","KIED"]);
      plt.(plotNm)('lgdNameList',obj.lgdNameList,'periodList',{[2003 2012]}); 
      close(gcf);
      plt.(plotNm)('periodList',{[1986 2000]}); 
      close(gcf); 
      clear plt 
            
      plt = Gplotmatrix(obj.workdirList,obj.model,"SourceFraction",["bio","ff","bb"]);
      plt.(plotNm)('lgdNameList',obj.lgdNameList,'periodList',{[2003 2012]}); 
      close(gcf);
      plt.(plotNm)('periodList',{[1986 2000]}); 
      close(gcf); 
      clear plt
      
      plt = Gplotmatrix(obj.workdirList,obj.model,...
        "Parameter",["life","KIEC","KIED"],"ECH4",["tot","anth_bio","anth_ff","bb"]);
      plt.cutUnit = []; 
      plt.(plotNm)('lgdNameList',obj.lgdNameList,'periodList',{[2003 2012]}); 
      close(gcf); 
      clear plt
    end

  end

  methods(Static)
    function addPriorUniform(model,categoryList)
      fig = gcf;  
      ax  = findobj(gcf,'type','ax');
      delete(legend)

      for i = 1:length(categoryList)
        name = categoryList(i);
        x = model.paraTargetList(1).(name).min + ...
          (model.paraTargetList(1).(name).max - ...
          model.paraTargetList(1).(name).min).*lhsdesign(1,100000);

        fig.CurrentAxes = ax(length(categoryList)-i+1); pause(0.3);

        h = histogram(x,'Normalization','Probability','DisplayStyle','stairs',...
          'LineWidth',2,'EdgeColor',[0.7 0.7 0.7]);

        h.BinEdges = model.paraTargetList(1).(name).min: ...
          (model.paraTargetList(1).(name).max-model.paraTargetList(1).(name).min)/20: ...
          model.paraTargetList(1).(name).max;

        ax(length(categoryList)-i+1).Children = [ax(length(categoryList)-i+1).Children(2:end);
        ax(length(categoryList)-i+1).Children(1)]; 

        if i == 1
          h = findobj(gca,'type','Histogram');
          inventoryList = reverse(extractBefore(reverse(model.Info.workdirList),"_")); 
          lgd = legend(flip(h),['Prior',strcat(inventoryList,"-Post")],'AutoUpdate','off');
          lgd.Position(1) = 0.35;
        end
      end
    end

  end
end

%%
  function y = plotCharacter(i,ax,fontSize,plotClass)
    switch i
      case 1,  txt = 'a';
      case 2,  txt = 'b';
      case 3,  txt = 'c';
      case 4,  txt = 'd';
      case 5,  txt = 'e';
      case 6,  txt = 'f';
      case 7,  txt = 'g';
      case 8,  txt = 'h';
      case 9,  txt = 'i';
      case 10, txt = 'j';
      case 11, txt = 'k';
      case 12, txt = 'l';
      case 13, txt = 'm';
      case 14, txt = 'n';
      case 15, txt = 'o';
      case 16, txt = 'p';
      case 17, txt = 'q';
      case 18, txt = 'r';
      case 19, txt = 's';
      case 20, txt = 't';
    end

    %switch plotClass
    %case 'Histogram', pos = [-0.15 1.2];
    %case 'Prctile',   pos = [0.025 0.92]; 
    %end
    pos = [0.025 0.92];

    y = text(ax,pos(1),pos(2),txt,'FontSize',fontSize,'FontWeight','bold',...
      'HorizontalAlignment','left','Units','normalized');
  end

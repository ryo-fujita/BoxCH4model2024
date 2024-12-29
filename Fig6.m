% Creates a plot of time series of prior and posterior CH4 emissions and 
% sinks by sector for 1850-2015 presented as Fig 6 in Fujita et al. 2025, JGR.
% Ryo Fujita 2024

classdef Fig6
  properties 
    workdirList
    model   
  end

  properties(Constant)
    statsNameList  = ["prc68","ave"];
  end

%%
  methods 
    function obj = Fig6(varargin)
      if nargin > 0
        obj.workdirList = varargin{1}; 
        obj.model = varargin{2}; 
      else
        obj.workdirList = "."; 
        obj.model = Model(".","Posterior",9);
      end
      plt = Prctile(obj.workdirList,obj.model,"ECH4",["tot","anth_bio","anth_ff"]);
      if length(obj.workdirList) >= 3, plt.setStatsNameList(Fig6.statsNameList); end

      fig_org = plt.plot('periodList',{[1850 2020]},'setSaveFig','no');
      Fig6.modifyTimeseries(fig_org,plt);
      fig = obj.addBBandLifeandNatural; pause(0.5)

      saveas(fig,'Fig6_BB_Life_Natr.png');
      savefig(fig,'Fig6_BB_Life_Natr.fig');
      close all
    end

%%
    function fig = addBBandLifeandNatural(obj)
      fig = gcf; pause(0.5); 
      ax = findobj(gcf,'type','ax');
      ax_add = obj.getAxBarPlotEms(["natr_bio","geo"],fig);     

      ax_add(1).Position(1) = ax(2).Position(1) + ax(2).Position(3) + 0.08;
      ax_add(1).Position(2) = ax(2).Position(2);
      ax_add(1).Position(3) = 0.15;
      ax_add(1).Position(4) = ax(2).Position(4);

      txt = findobj(ax_add(1).Children,'type','text'); delete(txt)
      tx  = text(ax_add(1),0.03,0.92,'e.','FontSize',12,'FontWeight', ...
        'bold','HorizontalAlignment','left','Units','normalized');
      tx  = text(ax_add(1),0.14,0.92,'Natural Biogenic','FontSize',12, ...
        'HorizontalAlignment','left','Units','normalized');
      ax_add(1).XTickLabel = ax_add(2).XTickLabel;

      ax_add(2).Position(1) = ax(2).Position(1) + ax(2).Position(3) + ax_add(1).Position(3) + 0.14;
      ax_add(2).Position(2) = ax(2).Position(2);
      ax_add(2).Position(3) = ax_add(1).Position(3);
      ax_add(2).Position(4) = ax(1).Position(4);

      txt = findobj(ax_add(2).Children,'type','text'); delete(txt)
      tx  = text(ax_add(2),0.03,0.92,'f.','FontSize',12,'FontWeight', ...
        'bold','HorizontalAlignment','left','Units','normalized');
      tx  = text(ax_add(2),0.14,0.92,'Geologic','FontSize',12, ...
        'HorizontalAlignment','left','Units','normalized');

      if length(ax_add(1).XTickLabel) >= 4 && length(ax_add(2).XTickLabel) >= 4
        for i = 1:2
          ax_add(i).XTickLabelRotation = 40;
          ax_add(i).XAxis.FontSize = 11;
          ax_add(i).XTickLabel = replace(ax_add(i).XTickLabel,'EDGAR','ED');
          ax_add(i).XRuler.TickLabelGapOffset = -1;
        end
      end

      figure(fig);
      an = findall(gcf,'type','textarrow');an.X = [0.045 0.045];

      plt = Prctile(obj.workdirList,obj.model,"ECH4life",["bb","ohAnomaly"]); 

      if length(obj.workdirList) >= 3
        plt.setStatsNameList(Fig6.statsNameList); 
      end

      fig_add = plt.plot('periodList',{[1850 2020]},'setSaveFig','no');
      delete(findobj(gcf,'type','legend'))

      ax1 = copyobj(fig_add.Children(2),fig); pause(0.5);
      ax2 = copyobj(fig_add.Children(1),fig); pause(0.5);

      figure(fig);
      fig.CurrentAxes = ax1; pause(0.5);      
      ax1.Position(1) = ax(3).Position(1) + ax(2).Position(3) + 0.08;
      ax1.Position(2) = ax(3).Position(2);
      ax1.Position(3) = 0.36;
      ax1.Position(4) = ax(3).Position(4);

      txt = findobj(ax1.Children,'type','text'); delete(txt)
      tx_d  = text(ax1,0.03,0.92,'d.','FontSize',12,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized');
      tx_d2 = text(ax1,0.11,0.92,'Biomass Burning','FontSize',12, ...
        'HorizontalAlignment','left','Units','normalized');

      ax1.XLabel.String = '';
      ax1.XTickLabel =  ax2.XTickLabel;
      ax1.YLim = [0 102];

      fig.CurrentAxes = ax2; pause(0.5)
      ax2.Position(1) = ax(1).Position(1) + ax(1).Position(3) + 0.08;
      ax2.Position(2) = ax(1).Position(2);
      ax2.Position(3) = 0.36;
      ax2.Position(4) = ax(1).Position(4);

      txt = findobj(ax2.Children,'type','text'); delete(txt)
      tx_f  = text(ax2,0.03,0.92,'g.','FontSize',12,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized');

      if ismember("life",plt.categoryList)
        txStr = "Total Lifetime";
        ax2.YLabel.String = 'Lifetime (yr)';
        ax2.YLabel.Position = [1.08 0.5 0];
        ax2.YLim = [6.9 10.7]; ax2.YTick = 7:0.5:11;
        ax2.YTickLabel = {'7','','8','','9','','10','','11'};
      elseif ismember("ohAnomaly",plt.categoryList) 
        ax2.YLabel.Position = [1.1 0.5 0];
        ax2.YLabel.String = 'OH anomaly (%)';
        txStr = "OH";
        ax2.YLim = [-22 22]; ax2.YTick = -20:10:20;
        ax2.YTickLabel = {'-20','-10','0','10','20'};
      end

      tx_f2 = text(ax2,0.11,0.92,txStr,'FontSize',12,'HorizontalAlignment','left', ...
        'Units','normalized');

      obj.addPriorLifetime(ax2,plt);
      fig.PaperPosition(3) = 27.5; fig.PaperSize(1) = fig.PaperPosition(3);
    end

%%
    function addPriorLifetime(obj,ax,plt)
      for i = 1:length(obj.workdirList)
        obj.model.workdir = obj.workdirList(i);
        obj.model  = obj.model.setup(obj.model.workdir);

        if ismember("life",plt.categoryList)
          h = plot(ax,[1850,2015.5], ...
            [obj.model.paraTargetList(1).life.def, ...
            obj.model.paraTargetList(obj.model.numNode-2).life.def],...
            ':','LineWidth',2,'Color',plt.modStyle.color(plt.lgdNameList(i),i));
        elseif ismember("ohAnomaly",plt.categoryList) 
          [yr_oh,ohAnomaly] = Loss.readOHAnomaly(obj.model); 
          h = plot(ax,yr_oh,ohAnomaly,':','LineWidth',2, ...
            'Color',plt.modStyle.color(plt.lgdNameList(i),i));
        end
      end
    end

%% 
    function ax = getAxBarPlotEms(obj,categoryList,fig)
      plt = BarPlot(obj.workdirList,obj.model,"ECH4",categoryList);
      plt.setSaveFig = 'no'; 

      fig_add = plt.plot('periodList',{[1850 2015]},'plotPriorStyle','single'); pause(0.5); 
      delete(legend)

      for j = length(fig_add.Children):-1:1
        i = length(fig_add.Children)-j+1;
        ax(i) = copyobj(fig_add.Children(j),fig); pause(0.5);

        ylabel(ax(i),[])
        if i ~= length(fig_add.Children), ax(i).XTickLabel = ''; end
        fig.CurrentAxes = ax(i); pause(0.5);

        switch categoryList(i)
        case 'natr_bio'
          ax(i).YLim(1) = 145; ax(i).YLim(2) = 535;
          ax(i).YTick = 150:50:550; ax(i).YTickLabel = {'','200','','300','','400','','500',''};
        case 'geo'
          ax(i).YLim(1) = -1; ax(i).YLim(2) = 95;
          ax(i).YTick = 0:10:90; ax(i).YTickLabel = {'0','','20','','40','','60','','80',''};
        case 'bb'
          ax(i).YLim(1) = -1; ax(i).YLim(2) = 72;
          ax(i).YTick = 0:10:70; ax(i).YTickLabel = {'0','','20','','40','','60',''};
        end
      end
    end
  end

  methods(Static)
    function fig = modifyTimeseries(fig,plt)
      figure(fig);
      % modify original text
      txt = findobj(gcf,'type','text');
      delete(txt);

      ax = findobj(gcf,'type','ax');
      ax(3).YLim = [180 750];
      ax(3).YTick = 200:50:700;
      ax(3).YTickLabel = {'200','','300','','400','','500','','600','','700',''};

      ax(2).YLim = [0 360];
      ax(2).YTick = 0:50:350;
      ax(2).YTickLabel = {'0','50','100','150','200','250','300','350'};

      ax(1).YLim  = [0 210];
      ax(1).YTick = 0:40:200;
      ax(1).YTickLabel = {'0','40','80','120','160','200'};

      tx_a  = text(ax(3),0.03,0.92,'a.','FontSize',12,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized');
      tx_a2 = text(ax(3),0.08,0.92,'Total','FontSize',12, ...
        'HorizontalAlignment','left','Units','normalized');

      tx_b  = text(ax(2),0.03,0.92,'b.','FontSize',12,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized');
      tx_b2 = text(ax(2),0.08,0.92,'Anthropogenic Biogenic','FontSize',12, ...
        'HorizontalAlignment','left','Units','normalized');

      tx_c  = text(ax(1),0.03,0.92,'c.','FontSize',12,'FontWeight','bold', ...
        'HorizontalAlignment','left','Units','normalized');
      tx_c2 = text(ax(1),0.08,0.92,'Anthropogenic Fossil Fuel','FontSize',12, ...
        'HorizontalAlignment','left','Units','normalized');

      an = findall(fig,'type','TextArrow');
      an.String = 'Global CH_{4} Emission (Tg CH_{4} yr^{-1})';
      an.Y = [0.68 0.68];

      lgd = findobj(fig,'Type','legend');

      if length(lgd.String) == 6 
        delete(lgd)
        for i = 1:length(plt.lgdNameList)
          text(ax(3),0.78,0.23-0.09*(i-1),plt.lgdNameList(i),'Units','normalized', ...
            'VerticalAlignment','bottom','Color',plt.modStyle.color(plt.lgdNameList(i)), ...
            'FontWeight','bold','FontSize',11.5);
        end
      else
        lgd.FontSize = 10;
      end

      % modify position
      ax(1).Position(3) = 0.5;
      ax(2).Position(3) = 0.5; ax(2).XLabel.String = '';
      ax(3).Position(3) = 0.5; ax(3).XLabel.String = '';

      ax(3).Position(2) = ax(3).Position(2) + 0.05;
      ax(2).Position(2) = ax(2).Position(2) + 0.025;
      ax(2).YAxisLocation = 'left';

      ax(1).Position(1) = 0.12; ax(1).Position(3) = 0.38;
      ax(2).Position(1) = 0.12; ax(2).Position(3) = 0.38;
      ax(3).Position(1) = 0.12; ax(3).Position(3) = 0.38;
      ax(1).Position(2) = ax(1).Position(2) -0.06;
      ax(2).Position(2) = ax(2).Position(2) -0.03;

      ax(2).XTickLabel = ax(1).XTickLabel;
      ax(3).XTickLabel = ax(1).XTickLabel;
    end
  end
end
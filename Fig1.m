% Creates a plot of timeseries of atmospheric CH4, d13C, dD, and D14C 
% presented as Fig 1 in Fujita et al. 2025.
% Ryo Fujita 2024

classdef Fig1
  properties 
  end
%%
  methods 
    function obj = Fig1(varargin)
      if nargin > 0
        workdirList = varargin{1}; 
        model = varargin{2}; 
      else
        workdirList = "."; 
        model = Model(".","Posterior",9);
      end
      lgdNameList = Plot.setlgdNameList(workdirList,model);
      periodList = [{[1740 2020]}];
      fig = gobjects(length(periodList),1);

      for i = 1:length(periodList)
        period = periodList(i);
        plt = Prctile(workdirList,model,"Atmosphere",["CH4","d13C","dD","D14C"]); 
        if length(workdirList) >= 3, plt.setStatsNameList(["prc68","ave"]); end 

        plt.plot('lgdNameList',lgdNameList,'periodList',period,'setSaveFig','no');
        fig(i) = Fig1.modifyFig1_final(period);

        figname = strcat('Fig1_',num2str(period{1}(1)),'-',num2str(period{1}(2)));
        saveas(fig(i),strcat(figname,'.png')); 
        savefig(fig(i),strcat(figname,'.fig'));
      end
    end

  end

  methods(Static)
    function placeFigIntoSubplot(iCol,axs)  
      ax   = findobj(axs,'type','ax'); pause(0.1);
      lgds = findobj(axs,'type','Legend'); pause(0.1)
      for i = 1:length(ax)
        ax(i).FontSize = 10;
        tx = findobj(ax(i),'Type','Text'); pause(0.1);
        for j = 1:length(tx), tx(j).FontSize = 10; end

        switch iCol
          case 1
            ax(i).Position(1) = 0.1;
            ax(i).Position(3) = 0.36;
            lgds(i).FontSize  = 5.5;
          case 2
            ax(i).Position(1) = 0.61;
            ax(i).Position(3) = 0.30;
            if i == 1, delete(lgds); end
            if ax(i).YLabel.Position(1) > 0, ax(i).YLabel.Position(1) = 1.06;  end
            if ax(i).YLabel.Position(1) < 0, ax(i).YLabel.Position(1) = -0.04; end
            if length(ax) == 4
              if i == 4
                fig = gcf; fig.CurrentAxes = ax(i);
                lgd = legend('Meinshausen17 (global)','NOAA (global)','FontSize',5.5,'Box','off');
                lgd.Position = [0.1234    0.8772    0.1577    0.0381];
              end

              switch tx.String
                case 'a', tx.String = 'e';
                case 'b', tx.String = 'f';
                case 'c', tx.String = 'g';
                case 'd', tx.String = 'h';
              end
            elseif  length(ax) == 3
              switch tx.String
                case 'a', tx.String = 'd';
                case 'b', tx.String = 'e';
                case 'c', tx.String = 'f';
              end
            end
        end
      end
    end  

%%
    function fig = modifyFig1_final(period)
      axes = findobj(gcf,'type','ax');

      switch length(axes)
      case 3
        tx_a = text(axes(3),0.03,0.93,'a','FontSize',12,'FontWeight','bold','HorizontalAlignment','left','Units','normalized');
        tx_b = text(axes(2),0.97,0.93,'b','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Units','normalized');
        tx_c = text(axes(1),0.03,0.93,'c','FontSize',12,'FontWeight','bold','HorizontalAlignment','left','Units','normalized');
      case 4
        tx_a = text(axes(4),0.03,0.93,'a','FontSize',12,'FontWeight','bold','HorizontalAlignment','left','Units','normalized');
        tx_b = text(axes(3),0.97,0.93,'b','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Units','normalized');
        tx_c = text(axes(2),0.03,0.93,'c','FontSize',12,'FontWeight','bold','HorizontalAlignment','left','Units','normalized');
        tx_d = text(axes(1),0.97,0.93,'d','FontSize',12,'FontWeight','bold','HorizontalAlignment','right','Units','normalized');
      end

      if period{1}(1) == 1740 && length(axes) == 4
        for i = 1:length(axes)
          ax = axes(i);
          switch ax.YLim(1)
          case 600
            ax.Position(2) = 0.65;
            ax.Position(4) = 0.19;
            
            lgd = findobj(gcf,'type','legend');
            if length(lgd) ~= 1
              delete(lgd(2))
              tx1 = text(ax,0.78,0.23,'CEDS','Units','normalized','VerticalAlignment','bottom','Color','#FF4B00','FontWeight','bold','FontSize',11.5);
              tx2 = text(ax,0.78,0.14,'EDGARv5','Units','normalized','VerticalAlignment','bottom','Color','#005AFF','FontWeight','bold','FontSize',11.5);
              tx3 = text(ax,0.78,0.05,'EDGARv6','Units','normalized','VerticalAlignment','bottom','Color','#03AF7A','FontWeight','bold','FontSize',11.5);
            end
          case -52.9
            ax.Position(2) = 0.46;
            ax.Position(4) = 0.19;
          case -130 
            ax.Position(2) = 0.27;
            ax.Position(4) = 0.19;
            ax.YLim = [-125 -65];
          case -300
            ax.Position(2) = 0.08;
            ax.Position(4) = 0.19;
            ax.YLim = [-280 480];
            for year = [1960 1970 1980 1990 2010 2020]
              tx = text(year, ax.YLim(1)-23, extractAfter(num2str(year),2),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',ax.XAxis.FontSize,'Color',ax.XAxis.Color);
            end
            fig = gcf;
            fig.Units = fig.PaperUnits;
            fig.Position = [fig.Position(1) fig.Position(2) fig.PaperPosition(3)  fig.PaperPosition(4)];
          end
        end
      end
    end

  end
end



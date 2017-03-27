%myzoom([0.26,0.2,0.6,0.3],[10,50,-5,15])

function myzoom(position,axiscale)
ha=get(gcf,'CurrentAxes');
ha1=copyobj(gca,gcf);
set(ha1,'position',position)
set(gcf,'CurrentAxes',ha1)
axis(axiscale)
xlabel('')
ylabel('')

set(gcf,'CurrentAxes',ha)
%set(gca,'LineWidth',1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')
set(gca,'LineWidth',1.5,'FontName','Arial','FontSize',10,'FontWeight','Bold')

clc;
clear;
load('315.mat')

figure
plot(log10(Cost1),'r-*','MarkerSize',5)
   
    set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A) ','FontName','Times New Roman','FontWeight','bold');
  export_fig 315.eps -painters -transparent 
   
clc;
clear;
x=[3,4,5,6,7,8,9,10,11,12,13,14,15 ]
y=[0.9603,5.2203,12.0448,20.4857,135.7004,448.3073, 970.873,1.101*10^4,4.796*10^4, 1.317*10^5, 1.719*10^6,   9.449 *10^6,  3.181*10^7]
    

plot(x,log10(y),'r-*','LineWidth', 1.5, 'MarkerSize',8)
grid on
 legend('m=3','FontName','Times New Roman','FontWeight','bold',0)
xlim( [ 3, 15] );
  set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('N','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A) ','FontName','Times New Roman','FontWeight','bold');
 legend
  export_fig 3.eps -painters -transparent 
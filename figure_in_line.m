clc;
clear;

x=[1,2,3,4,5,6,7,8];
y1=[8.7491,13.6576,20.5571,27.4025,34.298,40.9713,48.8001,54.6287]

y2=[499.99,660.2848,970.873,1324.6,1684,1980.9,2311,2657 ]

y3=[6.7661*10^4, 8.6388*10^4, 1.317*10^5 ,1.7313*10^5, 2.165*10^5,2.5917*10^5,3.0237*10^5,3.49*10^5 ]

y4=[1.7297*10^7,2.1168*10^7,3.181*10^7,4.247*10^7,5.5501*10^7,9.3771*10^7,2.43*10^8,3.1644*10^8]

plot(x,log10(y1),'r-*','MarkerSize',8,'LineWidth', 1.5)
hold on;
plot(x,log10(y2),'m:^','MarkerSize',8,'LineWidth', 1.5)
hold on;

plot(x,log10(y3),'k--o','MarkerSize',8,'LineWidth', 1.5)
hold on;
plot(x,log10(y4),'b-.x','MarkerSize',8,'LineWidth', 1.5)
hold on;

h=legend('N/M=2','N/M=3','N/M=4','N/M=5',0 );
set(h,'Fontsize',9);
   set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
   
    axis([0 9 0 10]);
%      myzoom([0.21,0.6 ,0.3,0.25],[1.9,2.1,13,55])
 % xlim( [ 1,10] );    
 xlabel('M','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A^*)','FontName','Times New Roman','FontWeight','bold');
  export_fig inte_line.eps -painters -transparent
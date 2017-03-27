clc;
clear;

x=[1,2,3,4,5]
y1=[0.927,18.951,993.562, 1.296*10^5,3.167*10^7]
   
y2=[0.96, 20.5571, 970.873,1.317*10^5,3.181*10^7]



a1=[1.228,25.151,1.294*10^3,1.716*10^5,4.22 *10^7 ]
a2=[1.2795,27.4025,1.3246*10^3,1.7313*10^5,4.247*10^7]

b1=[1.63, 31.466, 1.639*10^3,  2.143*10^5 ,5.381*10^7 ]
b2=[1.6852, 34.298, 1.684*10^3. 2.165*10^5,5.5501*10^7]


plot(x,log10(y1),'r-*','MarkerSize',8)
hold on;
plot(x,log10(y2),'c:^','MarkerSize',8)
hold on;

plot(x,log10(a1),'k-o','MarkerSize',8)
hold on;
plot(x,log10(a2),'b-.x','MarkerSize',8)
hold on;


plot(x,log10(b1),'g-s','MarkerSize',8)
plot(x,log10(b2),'m-.+','MarkerSize',8)
hold on;


% semilogy(x,y1,'r-*','MarkerSize',8,'LineWidth', 1.1)
% hold on;
% semilogy(x,y2,'c:^','MarkerSize',8,'LineWidth', 1.1)
% hold on;
% 
% semilogy(x,a1,'k-o','MarkerSize',8,'LineWidth', 1.1)
% hold on;
% semilogy(x,a2,'b-.x','MarkerSize',8,'LineWidth', 1.1)
% hold on;
% 
% 
% semilogy(x,b1,'g-s','MarkerSize',8,'LineWidth', 1.1)
% hold on;
% semilogy(x,b2,'m-.+','MarkerSize',8,'LineWidth', 1.1)
% hold on;

legend('M/\lambda (M=3)','NPGM (M=3)','M/\lambda (M=4)','NPGM (M=4)','M/\lambda (M=5)','NPGM (M=5)','Location','SouthEast' );
   set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
    myzoom([0.21,0.6 ,0.3,0.25],[1.9,2.1,1.1,1.7])
 axis([1 5 -0.05 8]);   
%     axis([1 5 0.47 10^8]);
%      myzoom([0.21,0.6 ,0.3,0.25],[1.9,2.1,13,55])
      
 xlabel('N/M','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A)_{min} ','FontName','Times New Roman','FontWeight','bold');
  export_fig inte.eps -painters -transparent 
   
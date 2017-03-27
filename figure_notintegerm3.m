
% 
%a=[ 3, 6, 9, 12, 15]
% b=[0.9603,24.3146,1174.48,1.515*10^5,3.693*10^7]

x=[3,4,5,6,7,8,9,10,11,12,13,14,15 ]
y=[0.9603,5.2203,12.0448,20.5571,135.7004,448.3073, 970.873,1.101*10^4,4.796*10^4, 1.317*10^5, 1.719*10^6,   9.449 *10^6,  3.181*10^7]
    

plot(x,log10(y),'r*','LineWidth', 1.5, 'MarkerSize',8)


x4=[4, 4];
y4=[0.9603,20.5571]     

x5= [5, 5,]
y5=[0.9603,20.5571];

x7=[7,7];
y7=[20.5571,970.873];

x8=[8,8];
y8=[20.55716,970.873];


x10=[10,10];
y10=[970.873,1.317*10^5,  ];


x11=[11,11];
y11=[970.873,1.317*10^5,  ];

x13=[13,13];
y13=[ 1.317*10^5, 3.181*10^7];

x14=[14,14];
y14=[ 1.317*10^5, 3.181*10^7];


% plot(a,log10(b),'b* ','MarkerSize',7)
% hold on

plot(x4,log10(y4),'b- ','LineWidth', 7)
hold on

%grid on

h1=plot(x5,log10(y5),'b- ','LineWidth', 7)
hold on

h2=plot(x7,log10(y7),'b- ','LineWidth', 7)
hold on

h3=plot(x8,log10(y8),'b- ','LineWidth', 7)
hold on

h4=plot(x10,log10(y10),'b- ','LineWidth', 7)
hold on

h5=plot(x11,log10(y11),'b- ','LineWidth', 7)
hold on

h6=plot(x13,log10(y13),'b- ','LineWidth', 7)
hold on

h7=plot(x14,log10(y14),'b- ','LineWidth', 7)
hold on


 set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
 
  set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
   set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
     set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
        set(get(get(h6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
           set(get(get(h7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% semilogy(a,(b),'bp ','MarkerSize',7)
% hold on
% 
% semilogy(x4,(y4),'b- ','LineWidth', 7)
% hold on
% 
% semilogy(x5,(y5),'b- ','LineWidth', 7)
% hold on
% 
% semilogy(x7,y7,'b- ','LineWidth', 7)
% hold on
% 
% semilogy(x8,(y8),'b- ','LineWidth', 7)
% hold on
% 
% semilogy(x10,(y10),'b- ','LineWidth', 7)
% hold on
% 
% semilogy(x11,(y11),'b- ','LineWidth', 7)
% hold on
% 
% semilogy(x13,(y13),'b- ','LineWidth', 7)
% hold on
% 
% semilogy(x14,(y14),'b- ','LineWidth', 7)
% hold on

%ylim( [ 10^(-1),  10^(-8)] );
plot(x,log10(y),'r*','LineWidth', 1.5, 'MarkerSize',12)


xlim( [ 3, 15] );
  
 xlabel('N','FontName','Times New Roman','FontWeight','bold');
 ylabel(' lg E(A)_{min} ','FontName','Times New Roman','FontWeight','bold');
 legend(['the estimated range of minimum cost  ',sprintf('\n'), 'when N/M is not  an integer'],['the minimum cost computed by NPGM'],'Location','SouthEast' )
 %,
 set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
  export_fig intenot.eps -painters -transparent 
  

%%simple gradient


clear all
 clc
format long 
%N=100

  
   M=2; N=6
%A=rand(N,N)

%   for k=1:N
%     A(k,k)=0
%   end
A=[zeros(1,N);eye(N-1),zeros(N-1,1)]; 
%A(1,N)=1

Topo=A;

A_diag=rand(N-1,1);

A=[zeros(1,N);diag(A_diag),zeros(N-1,1)]; 
%A(1,N)=rand(1,1);

  
A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);



%B=I(:,1:M)

 B=zeros(N,M);
 B(1,1)=1
 B(2,2)=1



I1=eye(M);


tf=1









v=0.002;


iterations=500
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  
Cost_cos=zeros(iterations,1);  



 
for ii=1:iterations
   
 
 H=zeros(N,N);

 ot=0.01; 
for k=1:tf/ot
    H=H+expm(-A*(ot*k))*B*B'*expm(-A'*(ot*k))*ot;
end
 C=pinv(H); 

   
Cost1(ii)=trace(C);

AA(:,:,ii)=A;
   
V1=eig(A);

V=esort(V1);

 II=eye(N);
for  iii=1:N-1
    
    if V(iii)==conj(V(iii+1)) & V(iii)~=V(iii+1)
        
        II(iii,iii)=0.5;
        
        II(iii,iii+1)=0.5;
           II(iii+1,iii)=0.5*i;
        
        II(iii+1,iii+1)=-0.5*i;
        
  
        
    end 
end




 chongshu=zeros(N,1);
ik=1;
while  ik>=1 && ik <=N
    index=find(V(ik:N)==V(ik));
    chongshu(ik)=length(index);
    ik=ik+length(index);
end 
 
 


BB=fliplr(vander(V));
CC1=exp(-V*tf);
for ik=1:N-1
    if chongshu(ik)>1
syms s
f=[];
for ppp=0:N-1
    f=[f,s^ppp];
end 
    

ff=f;
for pp=1:chongshu(ik)-1
    ff=[ff;diff(f,'s',pp)]; %%%求导
    CC1(ik+pp)=(-tf)^pp*exp(-V(ik)*tf);%%%求导
end 

fff=subs(ff,s,V(ik));      
BB(ik:ik+chongshu(ik)-1,:)=fff;%%%
 end 
end 



    

BBB=II*BB;
CCC=II*CC1;    


%alpha_tf=inv(BB)*CC1;
alpha_tf=BBB\CCC;

%alpha_tf =TVreg(BBB,CCC2,100)



F111=zeros(N,N);


for kk=1:tf/ot
    
   CC2=exp(-V*ot*kk); 
    
     for ik=1:N-1
    if chongshu(ik)>1


for pp=1:chongshu(ik)-1
     CC2(ik+pp)=(-ot*kk)^pp*exp(-V(ik)*(ot*kk));
end
    end
     end


CCC2=II*CC2;    

alpha_t=BBB\CCC2;


 pf=expm(-A*ot*kk); 

CB=C^2*pf*B*B';


F11=zeros(N,N);
for iiii =1:N-1
    
F1=zeros(N,N);
    for k=1:iiii
        
        F1=A'^(k-1)*CB*A'^(iiii-k)+F1;

    end

 F11=-2*alpha_t (iiii+1)*F1+F11;
end

F111=F11*ot+F111;


  
end 

    
    
  
      
 F=Topo.*F111/norm(Topo.*F111);
 
 F3=2*A;%(trace(A'*A)-N)*4*A; % gradient of norm function
 
 F3=F3/norm(F3);


F0=reshape((eye(N*N)-F3(:)*pinv(F3(:)))*F(:),N,N); 

Cost_cos(ii)=F(:)'*F3(:)/(norm(F(:))*norm(F3(:)));
A=A-v*F0/norm(F0);
%   F=F'+F-diag(diag(F));
%  for kkk=1:N
%     F(kkk,kkk)=0;
% end

%A=A-v*Topo.*F111/norm(F111);
A=(sqrt((N)/trace(A'*A)))*A; 


ii  

 
end




figure(1)
  plot((Cost1),'r-*')
  
  
  
figure(2)
  plot(Cost2,'r-*')
  
  
  figure(3)
plot(log10(Cost1),'r-*','MarkerSize',5)
   
    set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A) ','FontName','Times New Roman','FontWeight','bold');
  export_fig stem_converge.eps -painters -transparent 
  
    figure(4)
  plot(Cost_cos,'r-*','MarkerSize',5)
   
    set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A) ','FontName','Times New Roman','FontWeight','bold');
  export_fig stem_cos.eps -painters -transparent 
  
  
  weight=diag(A(2:N,1:N-1));
  weight=[weight;A(1,N)];
  figure(5)
  plot(weight)
  
    AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm(-A*tf);
AAAA-AL


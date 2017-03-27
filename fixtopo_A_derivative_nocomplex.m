clear all
 clc
format long 
%N=100

   M=1; N=2
   
%A=rand(N,N)

%   for k=1:N
%     A(k,k)=0
%   end
A=[zeros(1,N);eye(N-1),zeros(N-1,1)]; 

Topo=A;


NNN=rand(N-1,1);
for i=1:N-1
    A(i+1,i)=NNN(i);
end 
  
%A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);



B=I(:,1:M)

% B=zeros(N,M);
% B(1,1)=1
% B(3,2)=1



I1=eye(M);
tf=1

 
pf=expm(A*tf);


 

v=0.01;


iterations=1
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  



 
for ii=1:iterations
   
  
    pf=expm(A*tf);
    
 WB0=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
 C=pinv(WB0); 

   
Cost1(ii)=trace(C*pf*pf');

AA(:,:,ii)=A;

   
V=eig(A);


 II=eye(N);

for  iii=1:N-1
    
    if V(iii)==conj(V(iii+1)) & V(iii)~=V(iii+1)
        
        II(iii,iii)=0.5;
        
        II(iii,iii+1)=0.5;
           II(iii+1,iii)=0.5*i;
        
        II(iii+1,iii+1)=-0.5*i;
        
  
        
    end 
end

BB=fliplr(vander(V));




% if rank(diag(V)) < N;
% 
%  break
% else  
BBB=II*BB;

    
CC2=exp(V*tf);

CCC2=II*CC2;
%alpha_tf=pinv(BBB)*CCC2;
alpha_tf=BBB\CCC2;

%alpha_tf =TVreg(BBB,CCC2,100)


Cpf=C*pf;

F111=zeros(N,N);

ot=0.02;
for kk=1:tf/ot
    
CC1=exp(V*ot*kk);
CCC1=II*CC1;
%alpha_t=pinv(BBB)*CCC1;
alpha_t=BBB\CCC1;




CB=Cpf*pf'*C*expm(A*ot*kk)*B*B';


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

    
    
  F22=zeros(N,N);  
for iiii=1:N-1
  F2=zeros(N,N);
    for k=1:iiii
        
        F2=A'^(k-1)*Cpf*A'^(iiii-k)+F2;
    end

F22=2*alpha_tf(iiii+1)*F2+F22;
end
   
 F=F111+F22;   
      
% F=F/norm(F);
% 
% F3=(trace(A'*A)-N)*4*A; % gradient of norm function
% 
% F3=F3/norm(F3);
% 
% F0=reshape((eye(N*N)-F3(:)*pinv(F3(:)))*F(:),N,N); 


%   F=F'+F-diag(diag(F));
%  for kkk=1:N
%     F(kkk,kkk)=0;
% end

A=A-v*Topo.*F/norm(F);
A=(sqrt((N)/trace(A'*A)))*A; 


ii  

WB_1=zeros(N,N);
ot=0.02;
for k=1:tf/ot
     WB_1= WB_1+expm(A*(ot*k))*B(:,1)*B(:,1)'*expm(A'*(ot*k))*ot;
end

 
;
 
end


% for i=1:length(A)
%     for j=1:length(A)
%        if abs(A(i,j))<1.5*sum(sum(A))/(N*N)
%        A(i,j)=0;
%        end
%     
%     end
% end

figure(1)
  plot((Cost1),'r-*')
  
  
  
figure(2)
  plot(Cost2,'r-*')
  
  
  figure(3)
  plot(log10(Cost1),'r-*')
  
  
  AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm(A*tf);
AAAA-AL


%%%  分别求出B的每一列对应的特征值分布

%  WB0=zeros(N,N);
% ot=0.02;
% for k=1:tf/ot
%     WB0=WB0+expm(A*(ot*k))*B(:,1)*B(:,1)'*expm(A'*(ot*k))*ot;
% end
%  eig(WB0) 
%  


 WB_1=zeros(N,N);
ot=0.02;
for k=1:tf/ot
     WB_1= WB_1+expm(A*(ot*k))*B(:,1)*B(:,1)'*expm(A'*(ot*k))*ot;
end

 [X1 E1  Y1]=svd(WB_1)
 E1=diag(E1);
 
  E1(1)*X1(:,1)*X1(:,1)'+E1(2)*X1(:,2)*X1(:,2)'
 
  
   E1(1)*X1(:,1)*X1(:,1)'*E1(2)*X1(:,2)*X1(:,2)'
  
 
 WB_2=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    WB_2=WB_2+expm(A*(ot*k))*B(:,2)*B(:,2)'*expm(A'*(ot*k))*ot;
end
 [X2 E2 Y2]=svd(WB_2)
  E2=diag(E2);
  
 E2(1)*X2(:,1)*X2(:,1)'+E2(2)*X2(:,2)*X2(:,2)'
   
 
 
  acos(trace(X1(:,5)*X1(:,5)'*X2(:,1)*X2(:,1)'))*180/pi
   
  acos(trace(X1(:,1)*X1(:,1)'*X2(:,5)*X2(:,5)'))*180/pi
  
 
 acos(trace(X1(:,1)*X1(:,1)'*X1(:,2)*X1(:,2)'))*180/pi
   
 acos(trace(X1(:,1)*X1(:,1)'*X2(:,1)*X2(:,1)'))*180/pi
 
 acos(trace(X1(:,2)*X1(:,2)'*X2(:,2)*X2(:,2)'))*180/pi
 
 acos(trace(X1(:,2)*X1(:,2)'*X2(:,1)*X2(:,1)'))*180/pi
 
 acos(trace(X1(:,2)*X1(:,2)'*X2(:,2)*X2(:,2)'))*180/pi
  
 acos(trace(X2(:,1)*X2(:,1)'*X2(:,2)*X2(:,2)'))*180/pi
 
 
 
%   cost2(ii)=sum((eig(WB_1)-eig(WB_2)).^2)
 
 
 WB0=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
  [X3 E3 Y3]=svd(WB0) 
 
 
 
%%%  和一个点控制一条链的特征值分布基本相同
 A1=1.2217*[0 0 0 ; 1 0 0;  0 1 0];
 B1=[1 0 0]'
 
 
%   A1=[0 0 0 0; 1 0 0 0;  0 1 0 0;0 0 1 0];
%  B1=[1 0 0 0]'
  WB0=zeros(3,3);
ot=0.02;
for k=1:tf/ot
    WB0=WB0+expm(A1*(ot*k))*B1*B1'*expm(A1'*(ot*k))*ot;
end
 C=pinv(WB0); 
 
    pf=expm(A1*tf);

   
Cost0=trace(C*pf*pf');
 

  WW=zeros(3,3);
ot=0.02;
for k=1:tf/ot
    WW=WW+expm(-A1*(ot*k))*B1*B1'*expm(-A1'*(ot*k))*ot;
end
Cost00=trace(inv(WW));
 
 
 
 
 
 
 acos(trace(X1'*X2))*180/pi 
% acos(trace(X1'*X1))*180/pi 
% acos(trace(X1'*X3))*180/pi  
%   id = arrayfun(@(i)[i find(V == conj(V(i)))],[1:numel(V)]','UniformOutput',0);
%  id = cell2mat(id)



 
S1=zeros(N,N);
ot=0.02;
for k=1:tf/ot
     S1= S1+expm(-A*(ot*k))*B(:,1)*B(:,1)'*expm(-A'*(ot*k))*ot;
end
eig_S1=eig(S1)
 
 
 S2=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    S2=S2+expm(-A*(ot*k))*B(:,2)*B(:,2)'*expm(-A'*(ot*k))*ot;
end

eig_S2=eig(S2)

 S3=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    S3=S3+expm(-A*(ot*k))*B(:,3)*B(:,3)'*expm(-A'*(ot*k))*ot;
end

eig_S3=eig(S3)


 S0=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    S0=S0+expm(-A*(ot*k))*B*B'*expm(-A'*(ot*k))*ot;
end

eig_S0=eig(S0)


eig_SS=zeros(N,M);
for i=1:M
 SS=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    SS=SS+expm(-A*(ot*k))*B(:,i)*B(:,i)'*expm(-A'*(ot*k))*ot;
end

eig_SS(:,i)=eig(SS)

end





 D=diag(sum(A,1));
L=D-A;
eig(L)

eig(A'*A)

(sqrt(trace(inv(ds1)))*sqrt(trace(ds1))+sqrt(trace(inv(ds2)))*sqrt(trace(ds2))+sqrt(trace(inv(ds3)))*sqrt(trace(ds3)))^2/(trace(ds1)+trace(ds2)+trace(ds3))
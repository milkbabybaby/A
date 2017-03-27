clear all
 clc
format long 
%N=100

  
   M=1; N=4
%A=rand(N,N)

%   for k=1:N
%     A(k,k)=0
%   end
A=[zeros(1,N);eye(N-1),zeros(N-1,1)]; 

Topo=A;


  
A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);



B=I(:,1:M)

% B=zeros(N,M);
% B(1,1)=1
% B(3,2)=1



I1=eye(M);
tf=1

 
pf=expm(A*tf);


 

v=0.08;


iterations=200
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  



 
for ii=1:iterations
   
  
    pf=expm(A*tf);
    
 WB0=zeros(N,N);
ot=0.01;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
 C=pinv(WB0); 

   
Cost1(ii)=trace(C*pf*pf');

AA(:,:,ii)=A;
   
V1=eig(A);
V=sort(V1,'descend');

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
i=1;
while  i>=1 && i <=N
    index=find(V(i:N)==V(i));
    chongshu(i)=length(index);
    i=i+length(index);
end 
 
 


BB=fliplr(vander(V));





CC=exp(V*tf);
for i=1:N-1
    if chongshu(i)>1
syms s
f=[];
for ppp=0:N-1
    f=[f,s^ppp];
end 
    

ff=f;
for pp=1:chongshu(i)-1
    ff=[ff;diff(f,'s',pp)]; %%%求导
    CC(i+pp)=tf^pp*exp(V(i)*tf);%%%求导 tiaozheng 
end 

fff=subs(ff,s,V(i));      
BB(i:i+chongshu(i)-1,:)=fff;%%% tiaozheng
 end 
end 


BBB=II*BB;
CCC=II*CC;    



%alpha_tf=inv(BB)*CC1;
alpha_tf=BBB\CCC;

%alpha_tf =TVreg(BBB,CCC2,100)


Cpf=C*pf;

F111=zeros(N,N);

ot=0.01;
for kk=1:tf/ot
    
   CC2=exp(V*ot*kk); 
    
     for i=1:N-1
    if chongshu(i)>1


for pp=1:chongshu(i)-1
     CC2(i+pp)=(ot*kk)^pp*exp(V(i)*(ot*kk));
end
    end
     end

CCC2=II*CC2;  
alpha_t=BBB\CCC2;




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

 
end




figure(1)
  plot((Cost1),'r-*')
  
  
  
figure(2)
  plot(Cost2,'r-*')
  
  
  figure(3)
  plot(log10(Cost1),'r-*')
  
  weight=diag(A(2:N,1:N-1));
  figure(4)
  plot(weight)
  
    AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm(A*tf);
AAAA-AL


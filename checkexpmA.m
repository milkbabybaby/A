clear;
clc;
M=2; N=6
   

A=rand(N,N)
 
% for k=1:N
%     A(k,k)=0
% end
A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);


%B=rand(N,M); 
B=I(:,1:M)
% 
% B=[1 0 0 0 0 0; 0 0 0 0 0 1]';
%B=[0 1]';

 %B=[1 0; 0 0; 0 1; 0 0];
I1=eye(M);
tf=2

pf=expm(A*tf);
    
 WB0=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
 C=pinv(WB0); 


   
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


  AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm(A*tf);
AAAA-AL
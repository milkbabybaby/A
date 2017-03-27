clc;
clear;
N=4;
tf=1;
A=rand(N,N)
A=A*A';
V1=eig(A);

V=sort(V1,'descend');

 chongshu=zeros(N,1);
i=1;
while  i>=1 && i <=N
    index=find(V(i:N)==V(i));
    chongshu(i)=length(index);
    i=i+length(index);
end 
 
 


BB=fliplr(vander(V));
CC1=exp(sqrt(V)*tf);
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
    CC1(i+pp)=tf^pp*exp(sqrt(V(i))*tf);%%%求导
end 

fff=subs(ff,s,V(i));      
BB(i:i+chongshu(i)-1,:)=fff;%%%
 end 
end 



    



%alpha_tf=inv(BB)*CC1;
alpha_tf=BB\CC1;


    AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm((A^(1/2))*tf);
AAAA-AL
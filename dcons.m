function [dh,dg]=dcons(x,penal,ni,jds,dys)%K-Sλ��Լ����
dh=[];

[U,K,~,J,Umax]=prepare(x,penal,ni,jds,dys); 
load('qke.mat');    
x=(reshape(x,1,ni))';
duks=zeros(1,ni);
duks(1,1:ni)=0.;
dA(1,ni)=0.;

eh=1/Umax;
Aj(1:jds*3,1)=0.;
Aj(J,1)=1.;


for i=1:ni
   
        eval(['K1=qkeud_',num2str(i),';']); 
        eval(['K2=qkebalsa_',num2str(i),';']);
        dA(1:jds*3,i)=K\(penal*((x(i,1)^(penal-1))*(K1-K2))*U);
        eval(['clear K1',num2str(i),';']); 
        eval(['clear K1',num2str(i),';']);
   
end


for i=1:ni
     
     duks(1,i)=eh*Aj'*dA(1:jds*3,i);
end

dg3=duks;                   %%%%%���������ڿ��Ƶ���������

dg1(1:ni,1)=1.;
dg2(1:ni,1)=-1.;
dgg1=diag(dg1);
dgg2=diag(dg2);
dg=[dgg1;dgg2;dg3];
end
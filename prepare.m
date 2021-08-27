function [U,K,Uj,J,Umax]=prepare(x,penal,ni,jds,dys)       %求K-U

x=(reshape(x,1,ni))';
K=sparse(jds*3,jds*3);
qkeud=sparse(jds*3,jds*3);
qkebalsa=sparse(jds*3,jds*3);
load('bianjie.mat');
load('jiedian.mat');
load('danyuan.mat');
load('qke.mat');  
load('Nok.mat'); 
load('F.mat'); 
for i=1:ni
    eval(['qkeud=qkeud_',num2str(i),';']); 
    eval(['qkebalsa=qkebalsa_',num2str(i),';']);
    K=K+(qkebalsa+(x(i,1)^penal)*(qkeud-qkebalsa)); 
    eval(['clear qkeud_',num2str(i),';']); 
    eval(['clear qkebalsa_',num2str(i),';']);
end

K=K+Nok;
clear Nok;

for i=1:size(bianjie,1)
     K((bianjie(i,1)-1)*3+1,(bianjie(i,1)-1)*3+1)=K((bianjie(i,1)-1)*3+1,(bianjie(i,1)-1)*3+1)*(10^6);
     K((bianjie(i,1)-1)*3+2,(bianjie(i,1)-1)*3+2)=K((bianjie(i,1)-1)*3+2,(bianjie(i,1)-1)*3+2)*(10^6);
     K((bianjie(i,1)-1)*3+3,(bianjie(i,1)-1)*3+3)=K((bianjie(i,1)-1)*3+3,(bianjie(i,1)-1)*3+3)*(10^6);
end
K=full(K);
U=K\F;

 for i=1:jds*3
   U(i,1)=abs(U(i,1));
 end
 [Uj,J]=max(U);
Umax=7;        %最大位移



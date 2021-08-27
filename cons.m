function [h,g]=cons(x,penal,ni,jds,dys)  %K-S位移约束，设计变量区间约束
h=[];


[U,K,Uj,J,Umax]=prepare(x,penal,ni,jds,dys); 
x=(reshape(x,1,ni))';
rho=4.;

eh=exp(rho*(Uj/Umax));
uks=(1/rho)*log(eh);


g3=1-uks;
g1=x;
g2=1-x;
g4=[g1;g2];
g=[g4;g3];
end

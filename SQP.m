function [xt,mu,lam,val,k]=SQP(ni,penal,dys,jds,csz)
           rho=0.5; eta=0.1; sigma=0.3; epsilon1=1e-6;  %初始参数
%参数选取：
maxk=20.;
x0(1:ni,1)=csz;  

mu0=[];%等约束
lam0(1:(2*ni)+1,1)=0.;%不等约束


n=length(x0); l=length(mu0); m=length(lam0);
B0=eye(n);x=x0; mu=mu0;  lam=lam0; Bk=B0;

W = FEM(x,penal,ni,jds,dys) ;
dfk=df1(x,penal,ni,jds,dys);dfk=dfxz(dfk);
[hk,gk]=cons(x,penal,ni,jds,dys);
[Ae,Ai]=dcons(x,penal,ni,jds,dys); 
[Ae,Ai]=dfxc(Ae,Ai);
Ak=[Ae; Ai];

k=0;
while (k<maxk) 
    Wold=W;
    [dk,mu,lam]=SQPsub(dfk,Bk,Ae,hk,Ai,gk);%求解子问题
    deta=0.05; 
    tau=max(norm(mu,inf),norm(lam,inf));
    if(sigma*(tau+deta)<1)
        sigma=sigma;
    else
        sigma=1.0/(tau+2*deta);
    end
    im=0;   %%价值函数l
    while(im<=20)   % Armijo搜索
        if(phi1(x+rho^im*dk,sigma,penal,ni,jds,dys)-phi1(x,sigma,penal,ni,jds,dys)<eta*rho^im*dphi1(x,sigma,dk,penal,ni,jds,dys))
            mk=im;
            break;
        end
        im=im+1;
        if(im==20),    mk=10.;   end
    end
    alpha=rho^mk;
    x1=x+alpha*dk; 
    [hk,gk]=cons(x1,penal,ni,jds,dys); dfk=df1(x1,penal,ni,jds,dys);
    [Ae,Ai]=dcons(x1,penal,ni,jds,dys);  Ak=[Ae; Ai];
    lamu=pinv(Ak)'*dfk;  %计算最小二乘乘子
    if(l>0&&m>0)
        mu=lamu(1:l); lam=lamu(l+1:l+m);
    end
    if(l==0), mu=[]; lam=lamu;  end
    if(m==0), mu=lamu; lam=[];  end
    sk=alpha*dk;  %更新矩阵Bk
    yk=dlax(x1,mu,lam,penal,ni,jds,dys)-dlax(x,mu,lam,penal,ni,jds,dys);
    if(sk'*yk>0.2*sk'*Bk*sk)
        theta=1.0; 
    else
        theta=0.8*sk'*Bk*sk/(sk'*Bk*sk-sk'*yk);
    end
    zk=theta*yk+(1-theta)*Bk*sk;
    Bk=Bk+zk*zk'/(sk'*zk)-(Bk*sk)*(Bk*sk)'/(sk'*Bk*sk);
    x=x1; k=k+1;
   % xk=(reshape(x,1,ni))%%显示每次迭代x
    W = FEM(x,penal,ni,jds,dys);
    change = Wold-W;
    disp([' It.: ' sprintf('%4i',k) ' Obj.: ' sprintf('%10.4e',W) ...
        ' ch.: ' sprintf('%10.4e',change ) ])
  
    
    %收敛准则
    fz=abs(change);
    fm=Wold;
     if((fz/fm)<epsilon1)
         break; 
     end  
         
end
val=W;
xt=(reshape(x,1,ni))';
    




%%%%%%%% 价值函数 %%%%%%%
function p=phi1(x,sigma,penal,ni,jds,dys)
W = FEM(x,penal,ni,jds,dys) ;
f=W; [h,g]=cons(x,penal,ni,jds,dys); gn=max(-g,0);
l0=length(h);  m0=length(g);
if(l0==0), p=f+1.0/sigma*norm(gn,1); end
if(m0==0),  p=f+1.0/sigma*norm(h,1); end
if(l0>0&&m0>0)
    p=f+1.0/sigma*(norm(h,1)+norm(gn,1)); 
end

%%%%% 价值函数的方向导数%%%%%
function dp=dphi1(x,sigma,d,penal,ni,jds,dys)
df=df1(x,penal,ni,jds,dys);df=dfxz(df);
[h,g]=cons(x,penal,ni,jds,dys);  gn=max(-g,0);
l0=length(h);  m0=length(g);
if(l0==0),  dp=df'*d-1.0/sigma*norm(gn,1); end
if(m0==0), dp=df'*d-1.0/sigma*norm(h,1); end
if(l0>0&&m0>0)
        dp=df'*d-1.0/sigma*(norm(h,1)+norm(gn,1));
end

% %%%%%%%%% 拉格朗日函数 L(x,mu) %%%%%%%%%%%%%
 %%c = FEM(x,penal,ni,jds,dys);
 %f=c;
% [h,g]=cons(x,penal,ni,nj);  %调用约束函数文件
% l0=lemgth(h);  m0=length(g);
% if(l0==0), l=f-lam*g;  end
% if(m0==0),  l=f-mu'*h;  end
% if(l0>0&&m0>0)
 %        l=f-mu'*h-lam'*g;  
% end

%%%%%%%%% 拉格朗日函数的梯度 %%%%%%%%%%%%%
function dl=dlax(x,mu,lam,penal,ni,jds,dys)
df=df1(x,penal,ni,jds,dys);df=dfxz(df);
[Ae,Ai]=dcons(x,penal,ni,jds,dys);  %调用约束函数Jacobi矩阵文件
[Ae,Ai]=dfxc(Ae,Ai);
[m1,~]=size(Ai); [l1,~]=size(Ae);
if(l1==0),  dl=df-Ai'*lam;  end
if(m1==0), dl=df-Ae'*mu;  end
if(l1>0&&m1>0),  dl=df-Ae'*mu-Ai'*lam;  end

function dfk=dfxz(dfk)
function [Ae,Ai]=dfxc(Ae,Ai)









Qud=zeros(180,4);
Qbalsa=zeros(180,4);

load('q.mat');
load('xt.mat');

%%判断收敛
k=0.;
for i=1:180
    if(xt(i)>0.1&&xt(i)<0.9)
     k=k+1.;       
    else
     k=k+0.;
    end
end
%%%映射反演
value=0.25;
for i=1:180
    if(xt(i)>=value)
        xt(i)=1.;
    else
        xt(i)=0.;
    end
end
x=reshape(xt,6,30);

%建立两集合
for i=1:180
    if(xt(i)==1)
        Qud(i,1:4)=q(i,1:4);
        Qbalsa(i,1:4)=[0 0 0 0];
    else
        Qud(i,1:4)=[0 0 0 0];
        Qbalsa(i,1:4)=q(i,1:4);
    end
end
Qud(find(Qud==0))=[];
Qbalsa(find(Qbalsa==0))=[];
Qud=Qud';
Qbalsa=Qbalsa';
clear i q value xt
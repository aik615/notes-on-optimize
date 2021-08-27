%与ABAQUS信息进行验证

penal=1;
ni=180.;
jds=1323.;
load('xt.mat');
%load('x.mat');  xt=x;
%xt(1:ni,1)=1.;
%位移
x=(reshape(xt,1,ni))';
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

K=K+Nok;          %优化k+非优化k
clear  Nok;  %zlK ygK fbK bK;

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
 
%质量 
x=(reshape(x,1,ni))';
w=0.;
qmeud=sparse(24,24);
qmebalsa=sparse(24,24);
load('qme.mat');
load('Nom.mat');   
Now=Nom(1,1)*8.;
  for i=1:ni;
      eval(['qmeud=qmeud_',num2str(i),';']); 
      eval(['qmebalsa=qmebalsa_',num2str(i),';']);
      qmeud=qmeud(1,1)*8.;
      qmebalsa=qmebalsa(1,1)*8.;
      w=w+qmebalsa+x(i,1)^penal*(qmeud-qmebalsa);
      eval(['clear qmeud_',num2str(i),';']); 
      eval(['clear qmebalsa_',num2str(i),';']);
      
  end
    W=w+Now;
 clear bianjie danyuan F i jds jiedian K ni penal qkebalsa qkeud qmebalsa qmeud w x Nom Now U xt
 
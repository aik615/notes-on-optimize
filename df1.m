function df=df1(x,penal,ni,jds,dys)   
penal=1;
x=(reshape(x,1,ni))';
dc=zeros(ni,1);     
load('qme.mat');
  for i=1:ni;
      eval(['M1=qmeud_',num2str(i),';']); 
      eval(['M2=qmebalsa_',num2str(i),';']); 
      M1=M1(1,1)*800.;
      M2=M2(1,1)*800.;
      dc(i,1)=penal*(x(i,1)^(penal-1))*(M1-M2);
  end
df=reshape(dc,ni,1);
end


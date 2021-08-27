function W= FEM(x,penal,ni,jds,dys)
penal=1;
x=(reshape(x,1,ni))';
w=0.;
qmeud=sparse(24,24);
qmebalsa=sparse(24,24);
load('qme.mat');
load('Nom.mat');   
Now=Nom(1,1)*800.;
  for i=1:ni;
      eval(['qmeud=qmeud_',num2str(i),';']); 
      eval(['qmebalsa=qmebalsa_',num2str(i),';']);
      qmeud=qmeud(1,1)*800.;
      qmebalsa=qmebalsa(1,1)*800.;
      w=w+qmebalsa+x(i,1)^penal*(qmeud-qmebalsa);
      eval(['clear qmeud_',num2str(i),';']); 
      eval(['clear qmebalsa_',num2str(i),';']);
      
  end
    W=w+Now;
    
end





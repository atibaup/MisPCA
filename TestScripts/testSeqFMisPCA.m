close all
addpath(genpath('../Utilities'));

p=100;
n=30;
S=10;
dmax=30;
Spacing=10;
f_o=3;
SNR=30;
MissingTPratio=0;
r=p/10;

[X,SigmaAvOr,Fo,dini,Sigma,D]=GenerateFMisPCAData(p,n,S,f_o,SNR,r,dmax,Spacing,[],[]);

Masks=cell(S,1);
for i=1:S
   Masks{i}=ones(p,n);
end

[Options]=FMisPCA_options([],p,n,1);
Options.dmax=p-1;
Options.Nrndm=4;
Options.GridSize=Spacing;
Options.k=f_o;

[F,dopt,Lambdas,SigmaAv]=SeqFMisPCA(Sigma,Options);

[FPCA,Lambda]=eigs(SigmaAvOr,f_o);

[dtrue,Foracle,LambdasO]=OraclePCA(Sigma,f_o,dini);

figure;
subplot(3,1,1);
plot(F); hold on; 
title('FmisPCA Recovered Factors')
subplot(3,1,2);
plot(FPCA); hold on; 
title('PCA Recovered Factors')
subplot(3,1,3);
plot(Fo)
title('Original')

figure;
subplot(2,1,1)
imagesc(SigmaAvOr);
subplot(2,1,2)
imagesc(SigmaAv);

[Xhat1,A1,FactorsDistance1,MError1]=FMisPCA_Performance(X,D,Masks,Fo,F,dini,dopt,Lambdas,1);

[Xhat2,A2,FactorsDistance2,MError2]=FMisPCA_Performance(X,D,Masks,Fo,Foracle,dini,dtrue,LambdasO,1);

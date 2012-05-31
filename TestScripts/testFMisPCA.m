close all
addpath(genpath('../Utilities'));

p=100;
n=100;
S=10;
dmax=30;
Spacing=10;
f_o=1;
SNR=20;
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

% Joint PCA
[HPCA,LambdasPCA]=eigs(SigmaAvOr,f_o);
MinFactorsSubspaceDistance(Fo,HPCA)

% Oracle PCA
[dhatOr,FhatOr]=OraclePCA(Sigma,f_o,dini);
MinFactorsSubspaceDistance(Fo,FhatOr)

%             % MisPCA
F_AMisPCA=FMisPCA(Sigma,Masks,Options);
MinFactorsSubspaceDistance(Fo,F_AMisPCA)
%
%             % Seq PCA
F_SMisPCA=SeqFMisPCA(Sigma,Options);
MinFactorsSubspaceDistance(Fo,F_SMisPCA)

figure();
plot(HPCA,'-ob'); hold on;
plot(F_AMisPCA,'-xb'); hold on;
plot(F_SMisPCA,'-sg'); hold on;
plot(FhatOr,'-k')
legend('PCA','A MisPCA','Seq MisPCA','Oracle MisPCA')
            

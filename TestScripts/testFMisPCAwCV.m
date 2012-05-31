close all

addpath(genpath('../Utilities'));

p=50;
n=10;
S=10;
dmax=30;
Spacing=10;
f_o=2;
SNR=20;
MissingTPratio=0;
r=p/10;
        
[X,SigmaAvOr,Fo,dini,Sigma,D]=GenerateFMisPCAData(p,n,S,f_o,SNR,r,dmax,Spacing,[],[]);

[Options]=FMisPCA_options([],p,n,1);

Sigma=cell(S,1);
Masks=cell(S,1);
for i=1:S
   Sigma{i}=1/n* X{i}*X{i}';
   Masks{i}=ones(p,n);
end
Options.FiguresOnOff=0;
Options.GridSize=3;
Options.dmax=dmax;
Options.NCV=5;

[Results]=FMisPCAwCV(X,Options);

Genes=[];
for i=1:n
    Genes{i}=num2str(i);
end

[Err,h]=OriginalVSReconstructed(Results.d,Results.F,X,X,Masks,0,...
    [1:n],[1:S],1,1:S,'',Genes,zeros(S,1),Results.Lambdas,1:p,zeros(n,1));


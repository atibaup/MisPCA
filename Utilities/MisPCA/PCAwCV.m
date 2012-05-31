function [Results]=PCAwCV(X,Masks,Options)

S=length(X);
Sigma=cell(S,1);
for i=1:S
    Sigma{i}=cov(X{i}');
end

[XTr,SigmaTr,dummy,OmegaTt]=TrainingAndTest(X,Masks,Options);

BetaVec=logspace(-8,-3,Options.NpointsBeta);
tic;
CVErr=zeros(length(Options.kvec),Options.NpointsBeta);
for k=1:length(Options.kvec),
    Options.k=Options.kvec(k);
    for i=1:Options.NpointsBeta
        Options.beta=BetaVec(i);
        fprintf('CV: k=%d (%d/%d), i=%d/%d \n',Options.k,k,length(Options.kvec),i,Options.NpointsBeta)
        Err=zeros(Options.NCV,1);
        for r=1:Options.NCV
            [F,Lambdas]=PCA_w_Deflation(SigmaTr{r},Masks,Options);
            d=cell(k,1);
            for l=1:Options.k
               d{l}=zeros(Options.S,1); 
            end
            [Xhat]=RecoveredX(XTr{r},Masks,F,d,Lambdas);
            [dummy,Err(r)]=PredictionErr(X,Masks,OmegaTt{r},Xhat,0);
        end
        CVErr(k,i)=mean(Err);
    end
end
EllTime=toc;

[kopt,iopt]=find(CVErr==min(min(CVErr)));


Options.k=Options.kvec(kopt);
Options.beta=BetaVec(iopt);

fprintf('CV Results: f=%d (%d/%d), beta=%g (%d/%d) \n',Options.kvec(kopt),kopt,length(Options.kvec),Options.beta,iopt,Options.NpointsBeta)


d=cell(Options.k,1);
for l=1:Options.k
    d{l}=zeros(Options.S,1);
end
[F,Lambdas]=PCA_w_Deflation(Sigma,Masks,Options);

[Xhat,A,IndErr,RecErr]=RecoveredX(X,Masks,F,d,Lambdas);

Results.F=F;
Results.d=d;
Results.Options=Options;
Results.CVErr=CVErr;
Results.minCVErr=min(min(CVErr));
Results.kopt=kopt;
Results.iopt=iopt;
Results.beta=BetaVec(iopt);
Results.Xhat=Xhat;
Results.A=A;
Results.Lambdas=Lambdas;
Results.Err=RecErr;
Results.IndErr=IndErr;
Results.EllTime=EllTime;
if Options.FiguresOnOff
PlotMisPCAwCVResults(Results);
end
[Results.Xaligned,Results.Aaligned]=AlignedRepresentation(X,Masks,Results);
 

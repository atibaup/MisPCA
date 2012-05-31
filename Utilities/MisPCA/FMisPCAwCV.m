function [Results]=FMisPCAwCV(X,Options)


S=length(X);
n=size(X{1},2);
p=size(X{1},1);
Sigma=cell(S,1);
for i=1:S
    Sigma{i}=1/n*X{i}*X{i}';
end

if ~isfield(Options,'Masks')
    Masks=cell(S,1);
    for i=1:S
        Masks{i}=ones(p,n);
    end
else
    Masks=Options.Masks;
end

if ~isfield(Options,'beta')
    if Options.NpointsBeta>1
        BetaVec=logspace(-8,0,Options.NpointsBeta);
    else
        BetaVec=1e-3;
    end
else
    BetaVec=Options.beta;
    Options.NpointsBeta=length(BetaVec);
end

if (length(Options.kvec)==1 && Options.NpointsBeta==1)
    % No need to CValidate in this case:
    Options.NCV=1;
end
tic;
CVErr=zeros(length(Options.kvec),Options.NpointsBeta);
[XTr,SigmaTr,dummy,OmegaTt]=TrainingAndTest(X,Masks,Options);
fprintf('\n');
for k=1:length(Options.kvec)
    Options.k=Options.kvec(k);
    for i=1:Options.NpointsBeta
        Options.beta=BetaVec(i);
        fprintf('CV: k=%d (%d/%d), i=%d/%d \n',Options.k,k,length(Options.kvec),i,Options.NpointsBeta)

        Err=zeros(Options.NCV,1);
        if ~Options.FiguresOnOff && Options.NCV>1
            % Parallelizable for (deactivated for OCTAVE compatibility) 
 	    for r=1:Options.NCV
                [F,d,Lambdas]=FMisPCA(SigmaTr{r},Masks,Options);
                [Xhat]=RecoveredX(XTr{r},Masks,F,d,Lambdas);
                [dummy,Err(r)]=PredictionErr(X,Masks,OmegaTt{r},Xhat,0);
            end
        else
            for r=1:Options.NCV
                [F,d,Lambdas]=FMisPCA(SigmaTr{r},Masks,Options);          
                [Xhat]=RecoveredX(XTr{r},Masks,F,d,Lambdas);
                [dummy,Err(r)]=PredictionErr(X,Masks,OmegaTt{r},Xhat,0);
            end
            
        end
        CVErr(k,i)=mean(Err);
    end
end
fprintf('\n');
EllTime=toc;

[kopt,iopt]=find(CVErr==min(min(CVErr)));

Options.k=Options.kvec(kopt);
Options.beta=BetaVec(iopt);

fprintf('CV Results: f=%d (%d/%d), beta=%g (%d/%d) \n',Options.kvec(kopt),kopt,length(Options.kvec),Options.beta,iopt,Options.NpointsBeta)

[F,d,Lambdas]=FMisPCA(Sigma,Masks,Options);

[Xhat,A,IndErr,RecErr]=RecoveredX(X,Masks,F,d,Lambdas);

Results.F=F;
Results.d=d;
Results.Options=Options;
Results.CVErr=CVErr;
Results.minCVErr=min(min(CVErr));
Results.f=Options.kvec(kopt);
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

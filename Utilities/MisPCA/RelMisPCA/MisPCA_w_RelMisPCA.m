function [F,dopt,lambdamax,FadjAv]=MisPCA_w_RelMisPCA(Sigma,Masks,Options)
S=length(Sigma);
n=size(Masks{1},1);

k=Options.k;
beta=Options.beta;
dmax=Options.dmax;
FixedDelaysOnOff=Options.FixedDelaysOnOff;

if isfield(Options,'OP_OnOff')
    OP_OnOff=Options.OP_OnOff;
else
    OP_OnOff=0;
end

hini=1;
Niter=10;
tol=1e-3;
NrndmIni=Options.Nrndm;
Set=[-round(dmax/2)*ones(S,1),round(dmax/2)*ones(S,1)];

F=zeros(n,k);
lambdamax=zeros(1,k);
dopt=cell(k,1);

for i=1:k
    if FixedDelaysOnOff && i>1
        Set=[dopt{1},dopt{1}];
        dopt{i}=dopt{1};
        [dummy,F(:,i),dummy1,dummy2,lambdamax(i)]=RelaxationMisPCA_v7_fromSigma(Sigma,Masks,hini,Set,Niter,tol,NrndmIni,Options.FiguresOnOff,beta,Options.GridSize);
    else
        if i>1 && OP_OnOff
            Set=[dopt{1},Set(:,2)];
        end
        [dopt{i},F(:,i),dummy,dummy1,lambdamax(i)]=RelaxationMisPCA_v7_fromSigma(Sigma,Masks,hini,Set,Niter,tol,NrndmIni,Options.FiguresOnOff,beta,Options.GridSize);
    end
    [dummy,Sigma,FadjAv]=ComputeResidual_v3(dopt(1:i),F(:,1:i),lambdamax(1:i),Sigma,Masks);
end

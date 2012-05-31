function [F,dopt2,eigmax]=MisPCA_GivenDelays_fromCov(Sigma,Masks,dopt,k,beta,figOnOff)
S=length(Sigma);
n=size(Masks{1},1);


S=length(dopt);
if iscell(dopt)
    for i=1:S
        dopt2(i)=dopt{i};
    end
else
    dopt2=dopt;
end

SigmaAvIni=zeros(n);
opts.disp=0;
for i=1:S
    Indx=(Masks{i}(:,1)==1);
    SigmaAvIni(Indx,Indx)=SigmaAvIni(Indx,Indx)+1/S* Sigma{i}(Indx,Indx);
end

if figOnOff
    figure
    for j=1:S
        subplot(S,1,j)
        plot(1/max(max(Sigma{j}))*Sigma{j}); hold on;
    end
end
F=zeros(n,k);
for i=1:k
    %[ dopt{i} ] = EstimateDelays( Sigma,dmax );
    [ F(:,i),~, maxVal] = EstimateF_v2( dopt2,Sigma,Masks,beta );
%    [Res,Sigma]=ComputeResidual(d,F(:,1:i),X);
    %[Res]=ComputeResidual(d,F(:,1:i),X);
end
sigma_n=1/(n-1)*(trace(SigmaAvIni)-maxVal);
sigma_h=maxVal-sigma_n;
F=F;
eigmax=maxVal;
if figOnOff
    figure
    colors=hsv(k);
    for i=1:S
        subplot(S,1,i)
        plot(1/max(max(X{i}))*X{i}); hold on;
        Xr{i}=zeros(n,k);
        for j=1:k
            Fadj=circshift(F(:,j),dopt2(i));
            Xr{i}(:,j)=Fadj;
            plot(-1/max(max(abs(Fadj)))*Fadj,'Color',colors(j,:),'LineWidth',3);
        end
    end
end
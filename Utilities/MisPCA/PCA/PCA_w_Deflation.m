function [F,lambdamax,FadjAv]=PCA_w_Deflation(Sigma,Masks,Options)
S=length(Sigma);
n=size(Masks{1},1);


k=Options.k;
figOnOff=Options.FiguresOnOff;
beta=Options.beta;
dmax=Options.dmax;




F=zeros(n,k);

lambdamax=zeros(1,k);
dopt=cell(k,1);

for i=1:k
    [F(:,i),dopt{i},lambdamax(i)]=MisPCA_GivenDelays_fromCov(Sigma,Masks,zeros(S,1),1,beta,0);
    [Res,Sigma,FadjAv]=ComputeResidual_v3(dopt(1:i),F(:,1:i),lambdamax(1:i),Sigma,Masks);
    %[Sigma,SigmaAv,FadjAv]=ComputeResidual_v2(dopt(1:i),F(:,1:i),SigmaOr,lambdamax(1:i));
end
%dopt=dopt
% if figOnOff
%     figure
%     plot(lambdamax)
%     title('Estimated eigvals')
%     figure
%     colors=hsv(k);
%     for i=1:S
%         subplot(S,1,i)
%         plot(1/max(max(X{i}))*X{i}); hold on;
%         Xr{i}=zeros(n,k);
%         for j=1:k
%             dopt{j}
%             dopt{j}(i)
%             [T]=CreateTranslationMatrix(dopt{j}(i),n);
%             Xr{i}(:,j)=T*F(:,j);
%             plot(-1/max(max(abs(T*F(:,j))))*T*F(:,j),'Color',colors(j,:),'LineWidth',3);
%         end
%     end
% end
function [H,eta,SigmaAv,hest]=EstimateHnucNorm_v5(SigmaAv,Weights,beta)
%% Same as v4 only that here we only care about the elading eigenvector
p=size(SigmaAv,1);
%S=length(Sigma);
%Weights=zeros(p,1);
% SigmaAv=sparse(p,p);
% tic
% Niter=0;
% for i=1:S
%    for j= 1:GridSize:(Set(i,2)-Set(i,1))
%         if a((i-1)*p+j)>0
%             %Aux=sparse(p,p);
% %             Aux(Indx,Indx)=Sigma{i}(Indx,Indx);           
% %             d=(Set(i,1)+j-1);
% %             Aux=circshift(Aux',d)';
% %             Aux=circshift(Aux,d);            
%             %SigmaAv=SigmaAv+Aux;                    
%         end
%         Niter=Niter+1;
%     end
% end
% %looptime=toc
% InvWeights=1./Weights;

SigmaAv=diag(sqrt(Weights))*SigmaAv*diag(sqrt(Weights));

% s=[1 -1];
% v=sparse(1:length(s),ones(1,length(s)),s,p,1);
% v=1/norm(v)*v;
% W=circulant(v)';
e = ones(p,1);
W = spdiags([-e 2*e -e], -1:1, p, p);



B=SigmaAv-beta*(W'*W);

% %tic
[H,eta]=eigs(B,1,'lm');
% %toc

hest=sparse(H);
H=hest*hest';

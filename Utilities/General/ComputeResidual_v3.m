function [Res,Sigma,Fadj]=ComputeResidual_v3(dopt,F,Lambdas,Sigma,Masks)
S=length(Sigma);
n=size(Masks{1},1);
p=size(Sigma{1},2);
k=size(F,2);
plotindx=1;

%figure;

Res=cell(S,1);
for i=1:S
    Fadj=zeros(n,k);
    for j=1:k
         Fadj(:,j)= circshift(F(:,j),-dopt{j}(i));
    end
    warning off;
    Fadj=Fadj*diag(sqrt(Lambdas));
    
%     
%     subplot(S,2,plotindx)
%     [vec,eigvals]=eigs(Sigma{i},3);
%         plot(vec);
%         
%         
%     subplot(S,2,plotindx+1)
%     plotindx=plotindx+2;
        
        
    P=speye(size(Fadj,1))-Fadj*pinv(Fadj'*Fadj)*Fadj';

    Sigma{i}=P*Sigma{i}*P';
    Sigma{i}=1/2*(Sigma{i}+Sigma{i}');
%         [vec,eigvals]=eigs(Sigma{i},3);
%         plot(vec);

    
end

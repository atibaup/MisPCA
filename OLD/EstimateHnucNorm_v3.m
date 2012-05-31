function [H,C,SigmaAv]=EstimateHnucNorm_v3(a,Sigma,Masks,beta)
n=size(Masks{1},1);
S=length(Sigma);
SigmaAv=zeros(n);
for i=1:S
    for j=1:n
        if a((i-1)*n+j)>0
            T=CreateTranslationMatrix(j-1,n);
            Aux=zeros(n);
            Aux(Masks{i}==1,Masks{i}==1)=Sigma{i};
            Aux=T'*Aux*T;
%             TrsMask=T*Masks{i};
%             TrsMask=(TrsMask==1);
%             Aux=Aux(TrsMask,TrsMask);
            %TrsMask=(Masks{i}==1);
            %SigmaAv(TrsMask,TrsMask)=SigmaAv(TrsMask,TrsMask)+1/S*a((i-1)*n+j)*Aux;
            SigmaAv=SigmaAv+1/S*a((i-1)*n+j)*Aux;
        end
    end
end
cvx_quiet(true)
cvx_precision low;
if beta>0
    delta=1/2;
    W=eye(n)+delta*(diag(-ones(n-1,1),+1)+diag(-ones(n-1,1),-1));
    cvx_begin
        variable H(n,n) symmetric;
        minimize(TraceObjective_v2(H,Sigma,Masks,a)+beta*square_pos(norm(W*H,'fro')));%+beta*trace(H'*W'*W*H)
        subject to
             H== semidefinite(n);
            trace(H) == 1; %
            %norm(W*H,'fro')<=beta;
    cvx_end
else
    cvx_begin
        variable H(n,n) symmetric;
        minimize(TraceObjective_v2(H,Sigma,Masks,a));%+beta*trace(H'*W'*W*H)
        subject to
            H== semidefinite(n);
            trace(H) == 1; %
    cvx_end
end

C=-cvx_optval;
H=double(H);

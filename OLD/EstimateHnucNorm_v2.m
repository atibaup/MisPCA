function [H,C,SigmaAv]=EstimateHnucNorm_v2(a,Sigma,beta,GridSize,dmax)
n=size(Sigma{1},1);

Grid=0:GridSize:dmax;
nGrid=length(Grid);

S=length(Sigma);
SigmaAv=zeros(n);
for i=1:S
    for j=1:nGrid
        if a((i-1)*nGrid+j)>0
            T=CreateTranslationMatrix(Grid(j),n);
            SigmaAv=SigmaAv+1/S*a((i-1)*nGrid+j)*T'*Sigma{i}*T;
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
        minimize(TraceObjective(H,Sigma,a,GridSize,dmax));%+beta*trace(H'*W'*W*H)
        subject to
             H== semidefinite(n);
            trace(H) == 1; %
            norm(W*H,'fro')<=beta;
    cvx_end
else
    cvx_begin
        variable H(n,n) symmetric;
        minimize(TraceObjective(H,Sigma,a,GridSize,dmax));%+beta*trace(H'*W'*W*H)
        subject to
            H== semidefinite(n);
            trace(H) == 1; %
    cvx_end
end

C=-cvx_optval;
H=double(H);
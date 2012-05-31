function [H,C,SigmaAv]=EstimateHnucNorm(a,Sigma,beta)
n=size(Sigma{1},1);
S=length(Sigma);
C=0;
% H=sdpvar(n,n);
% W1=sdpvar(n,n,'symmetric');
% W2=sdpvar(n,n,'symmetric');
% R=sdpvar(2*n,2*n,'symmetric');
% SigmaAv=zeros(n);
% for i=1:S
%     for j=1:n
%         if a((i-1)*n+j)>0
%             T=CreateTranslationMatrix(j-1,n);
%             C=C-1/S*trace(H*T'*Sigma{i}*T);
%             SigmaAv=SigmaAv+1/S*T'*Sigma{i}*T;
%         end
%     end
% end
% C=C-beta*(trace(W1)+trace(W2));
% opts = sdpsettings('verbose',0);
% R=[W1,H;H',W2];
% F=set(H>0)+set(trace(H'*H) <=1)+set(R>0);
% sol=solvesdp(F,C,opts); %+beta*trace(H'*H)
% sol=sol
SigmaAv=zeros(n);
for i=1:S
    for j=1:n
        if a((i-1)*n+j)>0
            T=CreateTranslationMatrix(j-1,n);
            SigmaAv=SigmaAv+1/S*T'*Sigma{i}*T;
        end
    end
end

cvx_begin
    variable H(n,n) symmetric;
    variable W1(n,n) symmetric;
    variable W2(n,n) symmetric;
    expression R(2*n,2*n);
    R=[W1,H;H,W2];
    minimize(TraceObjective(H,Sigma,a)+beta*(trace(W1) + trace(W2)));
    subject to
        R== semidefinite(2*n);
        H== semidefinite(n);
        trace(H) == 1; % norm(H,'fro')<=1
cvx_end

C=-cvx_optval;
H=double(H);

function [tau,Sigma]=RelaxedMisPCA(X,tol,MaxIter)
S=length(X);
n=size(X{1},1);
p=size(X{1},2);

Sigma=cell(S,1);

W=1/sqrt(n)*dftmtx(n);
for i=1:S
    for j=i:S
        CenteredDatai=W*(X{i}-repmat(mean(X{i},2),1,p));
        CenteredDataj=W*(X{j}-repmat(mean(X{j},2),1,p));
        Sigma{i,j}=1/p*CenteredDatai*CenteredDataj';
        Sigma{j,i}=Sigma{i,j}';
    end
end

%% Projected gradient iterations
convergence=0;
tau=zeros(S,1);
t=1;
Cost=zeros(MaxIter,1);
alpha=.01;
while(~convergence)
    
    [Gradient ,Hessian]=RelaxedMisPCAGradient(tau, Sigma);
    tau=tau-Gradient
    tau=real(tau-Hessian\Gradient);
    CondNum=cond(Hessian)
    
    tau(tau<=0)=0;
    tau(tau>=2*pi)=2*pi;
    tau=tau;
    Cost(t)=0;
    for i=1:S
        for j=1:S
            U=diag(exp(-1i*2*pi/n*(0:n-1)*(tau(j)-tau(i))));
            if i~=j
                Cost(t)=Cost(t)+ real(trace(Sigma{i,j}*U));
            end
        end
    end
    if t>2
        if abs(Cost(t)-Cost(t-1))<=tol || t>=MaxIter
            convergence=1;
        end        
    end
    t=t+1;
end

figure
plot(Cost(1:t-1));

tau=tau/2*pi*n
function [Xr,X]=RecoverMatrixFromMissingData(X,MissingPerCnt)
S=length(X);
[n m]=size(X{1});
tol = 1e-8 ;
r=round(max(min(m,n)/10,1));
if MissingPerCnt>0
    Xr=cell(S,1);
    fprintf(1,'START: REcovering from missing data... \n') ;
    for i=1:S
        E = 1 - ceil( rand(n,m) - MissingPerCnt/sqrt(m*n)  ) ;
        M_E = sparse(X{i}.*E) ;
        [X1 S1 Y1 dist] = OptSpace(sparse(M_E),r,20,tol);
        Xr{i} = X1*S1*Y1';
    end
    fprintf(1,'End: Recovering from missing data \n') ;
else
    Xr=X;
end
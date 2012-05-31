function [V]=ComputeNullSpaceRankr(A,r)
[m,n] = size(A);
opts.tol=1e-8;
opts.issym=1;
opts.isreal=1;
warning off;
A=1/2*(A+A');
try
    [V,d]=eigs(A,r,'SA',opts);
catch
    opts.tol=1e-10;
    [V,d]=eigs(A,r,'SA',opts);
end
warning on;
d=diag(d);
tol = max(m,n) * max(d) * eps(class(A));
effr = find(d <= tol);

V=V(:,effr);
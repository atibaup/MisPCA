function [delays]=GenerateUniformDelays(dmax,spacing,S)

if round(S/(dmax+1))==S/(dmax+1)
    delays=kron(0:dmax,ones(1,S/(dmax+1)));    
else
    delays=-1;
    'error: S is not a multiple of dmax+1'
end
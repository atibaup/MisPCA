function [MinDist,Dist, PredDist,C_mtx,R_H,eigsRh]=SubspaceDistanceAndPredictedDistance(H,H_o,do)
F=size(H,2);
p=size(H,1);
MinDist=+inf;

if iscell(do)
do=cell2mat(do);
end

for i=1:p
    HHt=H*H';
    HHt=circshift(HHt,i);
    HHt=circshift(HHt',i);
    MinDist=min(MinDist,2*(F-trace((HHt)*(H_o*H_o'))));
end

Dist=2*(F-trace((H*H')*(H_o*H_o')));
dsorted=0:1:max(do);
dPoints=length(dsorted)

%compute r_h
r_h=zeros(F,F*dPoints);
r_h_=zeros(F,F*dPoints);

for i=1:dPoints
    r_h(:,(i-1)*F+1:i*F)=H_o'*circshift(H_o,dsorted(i));
    r_h_(:,(i-1)*F+1:i*F)=r_h(:,(i-1)*F+1:i*F)';
end


% generate block_toeplitz matrix
R_H=sparse(F*dPoints,F*dPoints);
for i=1:dPoints
    for j=i:dPoints
        Indx=(j-i)*F+1:((j-i+1)*F);
        Indxi=(i-1)*F+1:i*F;
        Indxj=(j-1)*F+1:j*F;
        R_H(Indxi,Indxj)=r_h(:,Indx);
        R_H(Indxj,Indxi)=r_h(:,Indx)';
    end
end

eigsRh=diag(D_RH);

C_mtx=V_B*D_RH_F*V_B';

PredDist=2*(F-trace(C_mtx));

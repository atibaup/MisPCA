addpath(genpath('../../Utilities'))
close all

%% f=1, flat signal
p=300;
n=1000;
SNR=40;
dmax=10;
Spacing=1;
H1=1/sqrt(p)*ones(p,1);
[X1,S,dummy,d1]=GenerateFMisPCADataUnifMisalignment(p,n,1,SNR,0,dmax,Spacing,[],H1);
[HPCA1,dummy]=eigs(S,1);
[MinDist1,Dist1, PredDist1,C_mtx1,R_H1]=SubspaceDistanceAndPredictedDistance(HPCA1,H1,d1);

%% f=3, orthogonal rectangular that stays orthogonal
W=3
H2=zeros(p,3);
H2(:,1)=[ones(W,1);zeros(p-W,1)];
H2(:,2)=circshift([ones(W,1);zeros(p-W,1)],31);
H2(:,3)=circshift([ones(W,1);zeros(p-W,1)],62);
H2=H2*diag(diag(1./sqrt(H2'*H2)));
dmax=50;
Spacing=2;
[X2,S,dummy,d2]=GenerateFMisPCADataUnifMisalignment(p,n,3,SNR,0,dmax,Spacing,[3,2,1],H2);
[HPCA2,dummy]=eigs(S,3);
[MinDist2,Dist2, PredDist2,C_mtx2,R_H2,eigsRH2]=SubspaceDistanceAndPredictedDistance(HPCA2,H2,d2);

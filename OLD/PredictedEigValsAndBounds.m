function [PredictedEval,eigMaxT, LB,UB]=PredictedEigValsAndBounds(h, theta, delays,c)
S=length(delays);
[T]=ConstructToeplitzfromDelays_NonRepeated(unique(delays),theta*h);
T=1/S*T;
[C]=CountsMatrix(delays);
ratio=S/(max(delays)+1);
rh=T(1,:);

%[PredictedEval,eigMaxT]=EigValTheoretical(diag(sqrt(C))*T*diag(sqrt(C)),c);
[PredictedEval,eigMaxT]=EigValTheoretical(S/(max(delays)+1)*T,c);

[LB,UB]=ToeplitzEigMaxBounds(rh);
LB=S/(max(delays)+1)*LB;
UB=S/(max(delays)+1)*UB;
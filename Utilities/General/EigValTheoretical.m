function [PredictedEval,eigMaxT]=EigValTheoretical(T,c)
eigMaxT=eigs(T,1,'lm'); %1/max(delays)*
if eigMaxT>sqrt(c)
    PredictedEval=(eigMaxT+1)*(1+c/eigMaxT);
else
    PredictedEval=(1+sqrt(c))^2;
end
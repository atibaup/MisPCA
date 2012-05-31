function [err,hestSh,hor]=ComparisonInterpSignals(hF,NoiseLessAligned,SamplingTimes,N)

h_est=interpft(real(hF),N);
hor=interpft(NoiseLessAligned(SamplingTimes),N);
corr1= xcorr(h_est,hor);
corr2= xcorr(-h_est,hor);
[maxVal1 indxMax1]=max(corr1);
[maxVal2 indxMax2]=max(corr2);
if maxVal1>maxVal2
indxMax=indxMax1-N;
sign=1;
else
    indxMax=indxMax2-N;
    sign=-1;
end
hestSh=sign*circshift(h_est,-indxMax);
err=abs(hestSh'*hor'/(norm(hestSh)*norm(hor)));
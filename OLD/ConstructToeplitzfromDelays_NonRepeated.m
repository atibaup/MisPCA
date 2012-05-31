function [T,hToeplitz]=ConstructToeplitzfromDelays_NonRepeated(delays,h)
S=length(delays);
p=length(h);
Maxd=max(abs(delays));
rh=xcorr(h);
rh=rh(p:p+Maxd);
hToeplitz=[];
for i=1:S
    hToeplitz=[hToeplitz, rh(delays(i)+1)];
end
T=sptoeplitz(sparse(hToeplitz),sparse(hToeplitz));
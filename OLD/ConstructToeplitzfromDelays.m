function [T,hToeplitz]=ConstructToeplitzfromDelays(delays,h)
S=length(delays);
p=length(h);
Maxd=max(abs(delays));
rh=xcorr(h);
rh=rh(p:p+Maxd);
hToeplitz=[];
for i=1:length(delays)
    hToeplitz=[hToeplitz, rh(delays(i))];
end
if length(delays)>0
    T=sptoeplitz(sparse(hToeplitz),sparse(hToeplitz));
else
   T=1; 
end
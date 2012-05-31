function [S]=CountsMatrix(delays)
delaysU=unique(delays);
dmax=length(delaysU);
S=zeros(dmax,1);
for i=1:dmax
    S(i)=length(find(delays==delaysU(i)));
end
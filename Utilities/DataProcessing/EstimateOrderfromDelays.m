function [Order,IntDelay]=EstimateOrderfromDelays(d,RefFactor)

S=length(d);
Avd=[];
for i=1:S
    Avd(i)=mean(d{i}(RefFactor));
end

[IntDelay, Order]=sort(Avd,'ascend');
function [Dist]=MinFactorsSubspaceDistance(Ho,H)
p=size(Ho,1);
F=size(Ho,2);
Dist=+inf;
for i=1:p
    Hc=circshift(H,i-1);
    Dist=min(Dist,2*(F-trace(Hc'*(Ho*Ho')*Hc)));
end

% figure; plot(Fo,'o'); hold on;
% plot(F,'x')
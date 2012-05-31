function [Dist]=FactorsSubspaceDistance(Fo,F)


P_Fo=Fo*((Fo'*Fo)\Fo');
P_F=F*((F'*F)\F');

Dist=norm(P_Fo-P_F,'fro');

% figure; plot(Fo,'o'); hold on;
% plot(F,'x')
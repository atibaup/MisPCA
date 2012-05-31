function [SnrTh,SnrTh1,SnrTh2]=FindThresholds(UB,LB,eigMaxT,SNR,c)

dUB=abs(mean(UB,2)-sqrt(c)).^2;
[minD1 indxmindUB]=min(dUB);
dLB=abs(mean(LB,2)-sqrt(c)).^2;
[minD2 indxmindLB]=min(dLB);
dreal=abs(mean(eigMaxT,2)-sqrt(c)).^2;
[minreal indxmindreal]=min(dreal);
SnrTh1=SNR(indxmindUB);
SnrTh2=SNR(indxmindLB);
SnrTh=SNR(indxmindreal);
% figure
% semilogy(SNR,eigMaxT,'-xk'); hold on;
% semilogy(SNR,UB,'-sr'); hold on;
% semilogy(SNR,LB,'-ob'); hold on;
% legend('eigMax','UB','LB')
% stem(SnrTh,log(max(eigMaxT)),'-xk');hold on;
% stem(SnrTh1,log(max(eigMaxT)),'-sr');hold on;
% stem(SnrTh2,log(max(eigMaxT)),'-ob');hold on;

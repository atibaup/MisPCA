
addpath(genpath('../../'));
close all;
clear all

p=20;
dmax=10;
K=7;
wVec=round(linspace(1,dmax-1,K));
h=figure; 

plotId=1;
for k=1:K
    s1=[ones(wVec(k),1);zeros(p-wVec(k),1)];
    s1=1/norm(s1)*s1;
    [lambda,lambda1]=ToeplitzBounds(s1,dmax);
    
    figure(h)
    subplot(K,2,plotId)
    plot(lambda,'-xb','LineWidth',2); hold on;
    plot(lambda1(:,1),'-^g'); hold on;
    plot(lambda1(:,2),'-vr'); hold on;
    title(['Rectangular, W=',num2str(wVec(k))],'FontSize',12)
    if k==1; leg=legend('True','UB','LB'); set(leg,'FontSize',12); end;
    axis tight;
    
    s2=[[wVec(k):-1:1]';zeros(p-wVec(k),1)];
    s2=1/norm(s2)*s2;
    
    
    [lambda,lambda1]=ToeplitzBounds(s2,dmax);
    
    subplot(K,2,plotId+1)
    plot(lambda,'-xb','LineWidth',2); hold on;
    plot(lambda1(:,1),'-^g'); hold on;
    plot(lambda1(:,2),'-vr'); hold on;
    axis tight;
   title(['Triangular, W=',num2str(wVec(k))],'FontSize',12)
    plotId=plotId+2;
end
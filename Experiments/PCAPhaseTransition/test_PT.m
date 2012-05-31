%close all
addpath(genpath('../../'));

F=3;
r=15; % pulse widht parameter (smaller than 15)
dmax=45; % Maximum delays
Spacing=10; % Size of delay grid.

Nrndm=30;

sigma=linspace(1,1,F);
sigmabar=1/max(sigma)*sigma;


test_PT_1;

test_PT_2;

test_PT_3;

%%

h=figure;
subplot(4,2,[1 2])
plot(H1,'LineWidth',2)
xlabel('Time')
ylabel('Magnitude')
leg=legend('Factor 1','Factor 2','Factor 3');
set(leg,'FontSize',12);

%%
colors=cool(F);
for f=1:F
    subplot(4,2,3)
    semilogy(SNR,mean(Lambda1(:,f,:),3),'--','color',colors(f,:),'LineWidth',2); hold on;
end

for f=1:F
    subplot(4,2,3)
%     for r=1:Nrndm
%         semilogy(SNR,Lambda1(:,f,r),'.','color',colors(f,:)); hold on;
%     end
    semilogy(SNR,PredictedLambda1(:,f),'-o','LineWidth',2,'color',colors(f,:)); hold on;
    xlabel('SNR (dB)','FontSize',12)
end

leg=legend('st eigvec','2nd eigvec','3rd eigvec');
set(leg,'FontSize',12);

for f=1:F    
    ylabel('\lambda_f(S(0))','FontSize',12)
    
    subplot(4,2,4)
    plot(SNR,mean(Corr1(:,f,:),3),'--','color',colors(f,:),'LineWidth',2); hold on;
%     for r=1:Nrndm
%         plot(SNR,Corr1(:,f,r),'.','color',colors(f,:)); hold on;
%     end
    plot(SNR,PredictedCorr1(:,f),'-o','LineWidth',2,'color',colors(f,:)); hold on;
    xlabel('SNR (dB)','FontSize',12)
    ylabel('|<h_f ,v_f(S(0)) >|^2','FontSize',12)
end
%

for f=1:F
    subplot(4,2,5)
    semilogy(SNR,mean(Lambda2(:,f,:),3),'--','color',colors(f,:),'LineWidth',2); hold on;
end

for f=1:F
    subplot(4,2,5)
%     for r=1:Nrndm
%         semilogy(SNR,Lambda2(:,f,r),'.','color',colors(f,:)); hold on;
%     end
    semilogy(SNR,PredictedLambda2(:,f),'-o','LineWidth',2,'color',colors(f,:)); hold on;
%    xlabel('SNR (dB)','FontSize',12)
    ylabel('\lambda_f(S(0))','FontSize',12)
    
    subplot(4,2,6)
    plot(SNR,mean(Corr2(:,f,:),3),'--','color',colors(f,:),'LineWidth',2); hold on;
%     for r=1:Nrndm
%         plot(SNR,Corr2(:,f,r),'.','color',colors(f,:)); hold on;
%     end
    plot(SNR,PredictedCorr2(:,f),'-o','LineWidth',2,'color',colors(f,:)); hold on;
    xlabel('SNR (dB)','FontSize',12)
    ylabel('|<h_f ,v_f(S(0)) >|^2','FontSize',12)
end



for f=1:F
    subplot(4,2,7)
    semilogy(SNR,mean(Lambda3(:,f,:),3),'--','color',colors(f,:),'LineWidth',2); hold on;
end


for f=1:F
    subplot(4,2,7)
%     for r=1:Nrndm
%         semilogy(SNR,Lambda3(:,f,r),'.','color',colors(f,:)); hold on;
%     end
    semilogy(SNR,PredictedLambda3(:,f),'-o','LineWidth',2,'color',colors(f,:)); hold on;
    xlabel('SNR (dB)','FontSize',12)
    ylabel('\lambda_f(S(0))','FontSize',12)
    
    subplot(4,2,8)
    plot(SNR,mean(Corr3(:,f,:),3),'--','color',colors(f,:),'LineWidth',2); hold on;
%     for r=1:Nrndm
%         plot(SNR,Corr3(:,f,r),'.','color',colors(f,:)); hold on;
%     end
    plot(SNR,PredictedCorr3(:,f),'-o','LineWidth',2,'color',colors(f,:)); hold on;
    xlabel('SNR (dB)','FontSize',12)
    ylabel('|<h_f ,v_f(S(0)) >|^2','FontSize',12)
end

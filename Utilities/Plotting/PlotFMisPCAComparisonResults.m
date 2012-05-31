function [h,h2]=PlotFMisPCAComparisonResults(Corr,CorrOrMisPCA,CorrMisPCA,CorrSeqMisPCA,SVec,F,SNR,Nrndm,XLabel)

NpointsS=length(SVec);

MarkerVec=['osx*v^'];
h=figure;
subaxis(2,3,1);
l=1;

plot(SNR,mean(CorrOrMisPCA(:,:,l),2),'-r','LineWidth',4); hold on;
plot(SNR,mean(CorrMisPCA(:,:,l),2),'-g','LineWidth',4); hold on;
plot(SNR,mean(Corr(:,:,l),2),'-b','LineWidth',4); hold on;
plot(SNR,mean(CorrSeqMisPCA(:,:,l),2),'-m','LineWidth',4); hold on;

if l==1
    legend('Oracle PCA','A-MisPCA','PCA','Seq-PCA')
end

plot(SNR,CorrMisPCA(:,:,l),'.g','MarkerSize',2); hold on;
plot(SNR,CorrOrMisPCA(:,:,l),'.r','MarkerSize',2); hold on;
plot(SNR,Corr(:,:,l),'.b','MarkerSize',2); hold on;
plot(SNR,CorrSeqMisPCA(:,:,l),'.m','LineWidth',4); hold on;
title([XLabel,'=',num2str(SVec(l))]);

subaxis(2,3,2);
l=round(length(SVec)/2);

plot(SNR,mean(CorrOrMisPCA(:,:,l),2),'-r','LineWidth',4); hold on;
plot(SNR,mean(CorrMisPCA(:,:,l),2),'-g','LineWidth',4); hold on;
plot(SNR,mean(Corr(:,:,l),2),'-b','LineWidth',4); hold on;
plot(SNR,mean(CorrSeqMisPCA(:,:,l),2),'-m','LineWidth',4); hold on;

plot(SNR,CorrMisPCA(:,:,l),'.g','MarkerSize',2); hold on;
plot(SNR,CorrOrMisPCA(:,:,l),'.r','MarkerSize',2); hold on;
plot(SNR,Corr(:,:,l),'.b','MarkerSize',2); hold on;
plot(SNR,CorrSeqMisPCA(:,:,l),'.m','LineWidth',4); hold on;

subaxis(2,3,3);
l=length(SVec);

plot(SNR,mean(CorrOrMisPCA(:,:,l),2),'-r','LineWidth',4); hold on;
plot(SNR,mean(CorrMisPCA(:,:,l),2),'-g','LineWidth',4); hold on;
plot(SNR,mean(Corr(:,:,l),2),'-b','LineWidth',4); hold on;
plot(SNR,mean(CorrSeqMisPCA(:,:,l),2),'-m','LineWidth',4); hold on;

plot(SNR,CorrMisPCA(:,:,l),'.g','MarkerSize',2); hold on;
plot(SNR,CorrOrMisPCA(:,:,l),'.r','MarkerSize',2); hold on;
plot(SNR,Corr(:,:,l),'.b','MarkerSize',2); hold on;
plot(SNR,CorrSeqMisPCA(:,:,l),'.m','LineWidth',4); hold on;

title([XLabel,'=',num2str(SVec(l))]);

% 
% 


ThresholdVec=2*F*linspace(.5,.1,3);
    
    h2=figure;
        
    h3=figure;
for t=1:length(ThresholdVec)
    %subplot(2,3,4+t-1);

    
    Threshold=ThresholdVec(t)

    SNRFineGrid=min(SNR):1:max(SNR);
    
    for l=1:NpointsS
         
        OrMisPCACurveAv=interp1(SNR,mean(Corr(:,:,l),2),SNRFineGrid);
        IndOrMisPCACurveAv=interp1(SNR,mean(CorrOrMisPCA(:,:,l),2),SNRFineGrid);
        MisPCACurve=interp1(SNR,mean(CorrMisPCA(:,:,l),2),SNRFineGrid);
        SeqMisPCACurveAv=interp1(SNR,mean(CorrSeqMisPCA(:,:,l),2),SNRFineGrid);

        IndxPCA=find(OrMisPCACurveAv<=Threshold,1);
        IndxOrMisPCA=find(IndOrMisPCACurveAv<=Threshold,1);
        IndxMisPCA=find(MisPCACurve<=Threshold,1);
        IndxSeqMisPCA=find(SeqMisPCACurveAv<=Threshold,1);
        
        if ~isempty(IndxPCA)
            SNRPCAAv2(l)=SNRFineGrid(IndxPCA(1));
        else
            SNRPCAAv2(l)=NaN;
        end
        if ~isempty(IndxOrMisPCA)
            SNROrMisPCAAv2(l)=SNRFineGrid(IndxOrMisPCA(1));
        else
            SNROrMisPCAAv2(l)=NaN;
        end
        if ~isempty(IndxMisPCA)
            SNRMisPCAAv2(l)=SNRFineGrid(IndxMisPCA(1));
        else
            SNRMisPCAAv2(l)=NaN;
        end
        
        if ~isempty(IndxSeqMisPCA)
            SNRSeqMisPCAAv2(l)=SNRFineGrid(IndxSeqMisPCA(1));
        else
            SNRSeqMisPCAAv2(l)=NaN;
        end
        figure(h3);
        subaxis(1,length(ThresholdVec),1+t-1);
        plot(SNRFineGrid,OrMisPCACurveAv,'r','Linewidth',2,'Marker',MarkerVec(t)); hold on;
        stem(SNRPCAAv2(l),1,'r','Linewidth',2,'Marker',MarkerVec(t));
        plot(SNRFineGrid,IndOrMisPCACurveAv,'g','Linewidth',2,'Marker',MarkerVec(t+1)); hold on;
        stem(SNROrMisPCAAv2(l),1,'g','Linewidth',2,'Marker',MarkerVec(t));
        plot(SNRFineGrid,MisPCACurve,'b','Linewidth',2,'Marker',MarkerVec(t+2)); hold on;
        stem(SNRMisPCAAv2(l),1,'b','Linewidth',2,'Marker',MarkerVec(t));
        plot(SNRFineGrid,SeqMisPCACurveAv,'m','Linewidth',2,'Marker',MarkerVec(t+3)); hold on;
        stem(SNRSeqMisPCAAv2(l),1,'m','Linewidth',2,'Marker',MarkerVec(t));
        grid on;

    end
    
 
    
    figure(h2);
    subaxis(1,length(ThresholdVec),1+t-1);
    plot(SVec,SNROrMisPCAAv2,'r','Linewidth',2,'Marker',MarkerVec(t)); hold on;
    plot(SVec,SNRMisPCAAv2,'g','Linewidth',2,'Marker',MarkerVec(t+1)); hold on;
    plot(SVec,SNRPCAAv2,'b','Linewidth',2,'Marker',MarkerVec(t+2)); hold on;
    plot(SVec,SNRSeqMisPCAAv2,'m','Linewidth',2,'Marker',MarkerVec(t+3)); hold on;
    grid on;
    if t==1
    legend('Oracle PCA','A-MisPCA','PCA','Seq-PCA')
    end


    axis([min(SVec),max(SVec),min(SNR),max(SNR)])
    title(['d(H,V_F(\Sigma)) <=',num2str(Threshold)]);
    
    xlabel(XLabel)
    if t==1; ylabel('SNR'); end;
end

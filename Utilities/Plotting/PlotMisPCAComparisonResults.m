function [h,h2]=PlotMisPCAComparisonResults(Corr,CorrOrMisPCA,CorrMisPCA,CorrSeqMisPCA,lambdaB,p,SVec,SNR,Nrndm,XLabel)

NpointsS=length(SVec);
NpointsSNR=length(SNR);

MarkerVec=['osx*v^'];
h=figure;
subplot(2,3,1);
l=1;

plot(SNR,mean(CorrOrMisPCA(:,:,l).^2,2),'-r','LineWidth',4); hold on;
plot(SNR,mean(CorrMisPCA(:,:,l).^2,2),'-g','LineWidth',4); hold on;
plot(SNR,mean(Corr(:,:,l).^2,2),'-b','LineWidth',4); hold on;
plot(SNR,mean(CorrSeqMisPCA(:,:,l).^2,2),'-m','LineWidth',4); hold on;

if l==1
    legend('Oracle PCA','A-MisPCA','PCA','Seq-PCA')
end

plot(SNR,CorrMisPCA(:,:,l).^2,'.g','MarkerSize',2); hold on;
plot(SNR,CorrOrMisPCA(:,:,l).^2,'.r','MarkerSize',2); hold on;
plot(SNR,Corr(:,:,l).^2,'.b','MarkerSize',2); hold on;
plot(SNR,CorrSeqMisPCA(:,:,l).^2,'.m','LineWidth',4); hold on;
title([XLabel,'=',num2str(SVec(l))]);

subplot(2,3,2);
l=round(length(SVec)/2);

plot(SNR,mean(CorrOrMisPCA(:,:,l).^2,2),'-r','LineWidth',4); hold on;
plot(SNR,mean(CorrMisPCA(:,:,l).^2,2),'-g','LineWidth',4); hold on;
plot(SNR,mean(Corr(:,:,l).^2,2),'-b','LineWidth',4); hold on;
plot(SNR,mean(CorrSeqMisPCA(:,:,l).^2,2),'-m','LineWidth',4); hold on;

plot(SNR,CorrMisPCA(:,:,l).^2,'.g','MarkerSize',2); hold on;
plot(SNR,CorrOrMisPCA(:,:,l).^2,'.r','MarkerSize',2); hold on;
plot(SNR,Corr(:,:,l).^2,'.b','MarkerSize',2); hold on;
plot(SNR,CorrSeqMisPCA(:,:,l).^2,'.m','LineWidth',4); hold on;

subplot(2,3,3);
l=length(SVec);

plot(SNR,mean(CorrOrMisPCA(:,:,l).^2,2),'-r','LineWidth',4); hold on;
plot(SNR,mean(CorrMisPCA(:,:,l).^2,2),'-g','LineWidth',4); hold on;
plot(SNR,mean(Corr(:,:,l).^2,2),'-b','LineWidth',4); hold on;
plot(SNR,mean(CorrSeqMisPCA(:,:,l).^2,2),'-m','LineWidth',4); hold on;

plot(SNR,CorrMisPCA(:,:,l).^2,'.g','MarkerSize',2); hold on;
plot(SNR,CorrOrMisPCA(:,:,l).^2,'.r','MarkerSize',2); hold on;
plot(SNR,Corr(:,:,l).^2,'.b','MarkerSize',2); hold on;
plot(SNR,CorrSeqMisPCA(:,:,l).^2,'.m','LineWidth',4); hold on;

title([XLabel,'=',num2str(SVec(l))]);

% 
% 

h2=figure;
ThresholdVec=linspace(.15,.8,3);

for t=1:length(ThresholdVec)
    %subplot(2,3,4+t-1);
    subplot(1,length(ThresholdVec),1+t-1);
    
    Threshold=ThresholdVec(t)

    SNRPCA=zeros(NpointsS,Nrndm);
    SNROrMisPCA=zeros(NpointsS,Nrndm);
    SNRMisPCA=zeros(NpointsS,Nrndm);
    SNRPredPCA=zeros(NpointsS,1);
    SNRPredOracle=zeros(NpointsS,1);
    SNRFineGrid=min(SNR):.5:max(SNR);
    
    for l=1:NpointsS
        [SNRPredPCA(l)]=ComputePredictedSNRPoint(sqrt(p/SVec(l)),Threshold,lambdaB(l));
        [SNRPredOracle(l)]=ComputePredictedSNRPoint(sqrt(p/SVec(l)),Threshold,1);
        
        OrMisPCACurveAv=zeros(length(SNRFineGrid),1);
        MisPCACurveAv=zeros(length(SNRFineGrid),1);
        SeqMisPCACurveAv=zeros(length(SNRFineGrid),1);
        IndOrMisPCACurveAv=zeros(length(SNRFineGrid),1);
        
        for i=1:Nrndm
                       
            %OrMisPCACurve=interp1(SNR,mean(Corr(:,:,l).^2,2),SNRFineGrid);
            OrMisPCACurve=interp1(SNR,Corr(:,i,l).^2,SNRFineGrid);

            OrMisPCACurveAv=OrMisPCACurveAv+1/Nrndm*OrMisPCACurve';
            IndxPCA=find(OrMisPCACurve>=Threshold);
            
            if isempty(IndxPCA); IndxPCA=length(SNRFineGrid);end;
            
            %IndOrMisPCACurve=interp1(SNR,mean(CorrOrMisPCA(:,:,l).^2,2),SNRFineGrid);
            IndOrMisPCACurve=interp1(SNR,CorrOrMisPCA(:,i,l).^2,SNRFineGrid);
            IndOrMisPCACurveAv=IndOrMisPCACurveAv+1/Nrndm*IndOrMisPCACurve';
            IndxOrMisPCA=find(IndOrMisPCACurve>=Threshold);
            
            if isempty(IndxOrMisPCA); IndxOrMisPCA=length(SNRFineGrid);end;
            
            %MisPCACurve=interp1(SNR,mean(CorrMisPCA(:,:,l).^2,2),SNRFineGrid);
            MisPCACurve=interp1(SNR,CorrMisPCA(:,i,l).^2,SNRFineGrid);
            MisPCACurveAv=MisPCACurveAv+1/Nrndm*MisPCACurve';
            IndxMisPCA=find(MisPCACurve>=Threshold);

            
            if isempty(IndxMisPCA); IndxMisPCA=length(SNRFineGrid);end;
            
            SeqMisPCACurve=interp1(SNR,CorrSeqMisPCA(:,i,l).^2,SNRFineGrid);
            SeqMisPCACurveAv=SeqMisPCACurveAv+1/Nrndm*SeqMisPCACurve';
            IndxSeqMisPCA=find(SeqMisPCACurve>=Threshold);

            
            if isempty(IndxSeqMisPCA); IndxSeqMisPCA=length(SNRFineGrid);end;
            
            SNRPCA(l,i)=SNRFineGrid(IndxPCA(1));
            SNROrMisPCA(l,i)=SNRFineGrid(IndxOrMisPCA(1));
            SNRMisPCA(l,i)=SNRFineGrid(IndxMisPCA(1));
            SNRSeqMisPCA(l,i)=SNRFineGrid(IndxSeqMisPCA(1));
        end
        
        OrMisPCACurveAv=interp1(SNR,mean(Corr(:,:,l).^2,2),SNRFineGrid);
        IndOrMisPCACurveAv=interp1(SNR,mean(CorrOrMisPCA(:,:,l).^2,2),SNRFineGrid);
        MisPCACurve=interp1(SNR,mean(CorrMisPCA(:,:,l).^2,2),SNRFineGrid);
        SeqMisPCACurveAv=interp1(SNR,mean(CorrSeqMisPCA(:,:,l).^2,2),SNRFineGrid);
        
%         figure;
%         plot(OrMisPCACurveAv,'-b','LineWidth',2);hold on;
%         plot(IndOrMisPCACurveAv,'-r','LineWidth',2);hold on;
%         plot(MisPCACurve,'-g','LineWidth',2);hold on;
%         plot(SeqMisPCACurveAv,'-m','LineWidth',2);hold on;
%         legend('PCA','Oracle','A-MisPCA','Seq-MisPCA')
        

        IndxPCA=find(OrMisPCACurveAv>=Threshold,1);
        IndxOrMisPCA=find(IndOrMisPCACurveAv>=Threshold,1);
        IndxMisPCA=find(MisPCACurve>=Threshold,1);
        IndxSeqMisPCA=find(SeqMisPCACurveAv>=Threshold,1);
        
        if ~isempty(IndxPCA)
            SNRPCAAv2(l)=SNRFineGrid(IndxPCA);
        else
            SNRPCAAv2(l)=NaN;
        end
        if ~isempty(IndxOrMisPCA)
            SNROrMisPCAAv2(l)=SNRFineGrid(IndxOrMisPCA);
        else
            SNROrMisPCAAv2(l)=NaN;
        end
        if ~isempty(IndxMisPCA)
            SNRMisPCAAv2(l)=SNRFineGrid(IndxMisPCA);
        else
            SNRMisPCAAv2(l)=NaN;
        end
        
        if ~isempty(IndxSeqMisPCA)
            SNRSeqMisPCAAv2(l)=SNRFineGrid(IndxSeqMisPCA);
        else
            SNRSeqMisPCAAv2(l)=NaN;
        end

    end
    
    
    
    SNRPCAAv=median(SNRPCA,2);
    SNROrMisPCAAv=median(SNROrMisPCA,2);
    SNRMisPCAAv=median(SNRMisPCA,2);
    SNRSeqMisPCAAv=median(SNRSeqMisPCA,2);

    plot(SVec,SNROrMisPCAAv2,'r','Linewidth',2,'Marker',MarkerVec(t)); hold on;
    plot(SVec,SNRMisPCAAv2,'g','Linewidth',2,'Marker',MarkerVec(t+1)); hold on;
    plot(SVec,SNRPCAAv2,'b','Linewidth',2,'Marker',MarkerVec(t+2)); hold on;
    plot(SVec,SNRSeqMisPCAAv2,'m','Linewidth',2,'Marker',MarkerVec(t+3)); hold on;
    
    if t==1
    legend('Oracle PCA','A-MisPCA','PCA','Seq-PCA')
    end
        
    %plot(SVec,SNRPredOracle,'-.r','Linewidth',2,'Marker','*'); hold on;
%     plot(SVec,SNROrMisPCA,'.r'); hold on;
%     plot(SVec,SNRMisPCA,'.g'); hold on;
%     plot(SVec,SNRPCA,'.b'); hold on;
%     plot(SVec,SNRSeqMisPCA,'.m'); hold on;
    
    %plot(SVec,SNRPredPCA,'-.b','Linewidth',2,'Marker','*'); hold on;

    %if ~isempty(SNRstarPred)>0

    %end

    axis([min(SVec),max(SVec),min(SNR),max(SNR)])
    title(['<. , h> =',num2str(Threshold)]);
    
    xlabel(XLabel)
    ylabel('SNR')
end

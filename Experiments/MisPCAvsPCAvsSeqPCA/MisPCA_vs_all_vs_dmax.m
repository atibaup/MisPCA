clear all;
close all;

addpath(genpath('../../Utilities'))

NCPU=6;
if matlabpool('size') ~= NCPU
    if  matlabpool('size')>0
        matlabpool close;
    end
    matlabpool('open', NCPU)
end

Nrndm=NCPU*4;
p=100; % time dimension
n=10; % samples per replicate
F=2;

NpointsSNR=20;
NpointsS=10;

SNR=round(linspace(-10,30,NpointsSNR));
r=round(p/30); 
S=30; 

dmaxVec=round(linspace(1,p-1,NpointsS));
dmaxVec=sort(dmaxVec,'ascend');

Masks=cell(S,1);
for i=1:S
    Masks{i}=ones(p,n);
end

[dummy1,dummy2,H]=GenerateFMisPCAData(p,n,1,F,SNR(1),r,1,1,[],[]);

Corr=zeros(NpointsSNR,Nrndm,NpointsS);
CorrMisPCA=zeros(NpointsSNR,Nrndm,NpointsS);
CorrOrPCA=zeros(NpointsSNR,Nrndm,NpointsS);
CorrSeqPCA=zeros(NpointsSNR,Nrndm,NpointsS);
TimePCA=zeros(NpointsSNR,Nrndm,NpointsS);
TimeMisPCA=zeros(NpointsSNR,Nrndm,NpointsS);
TimeOrPCA=zeros(NpointsSNR,Nrndm,NpointsS);
TimeSeqPCA=zeros(NpointsSNR,Nrndm,NpointsS);


for l=1:NpointsS
    dmax=dmaxVec(l);
    fprintf('dmax=%d (%d/%d)\n',dmax,l,NpointsS)
    [Options]=Simulations_options(p,max(ceil(dmax/5),1),dmax);
    Options.k=F;
    Spacing=Options.GridSize;    
    for k=1:NpointsSNR
        
        fprintf('\t SNR=%g (%d/%d)\n',SNR(k),k,NpointsSNR)
        
        % Parallelizable for (deactivated for OCTAVE compatibility) 
 for i=1:Nrndm
            % Generate misaligned signals with random delays
            [X,SigmaAv,Fo,dini,Sigma,D]=GenerateFMisPCAData(p,n,S,F,SNR(k),r,dmax,Spacing,[],H);

            % Joint PCA
            tic;
            [HPCA,LambdasPCA]=eigs(SigmaAv,F);   
            TimePCA(k,i,l)=toc;
            Corr(k,i,l)=MinFactorsSubspaceDistance(Fo,HPCA);

            
            % Oracle PCA
            tic;
            [dhatOr,FhatOr]=OraclePCA(Sigma,F,dini);
            TimeOrPCA(k,i,l)=toc;
            CorrOr(k,i,l)=MinFactorsSubspaceDistance(Fo,FhatOr);

%             % MisPCA
            tic;
            F_AMisPCA=FMisPCA(Sigma,Masks,Options);         
            TimeMisPCA(k,i,l)=toc;
            CorrMisPCA(k,i,l)=MinFactorsSubspaceDistance(Fo,F_AMisPCA);

%             % Seq PCA
            tic;
            F_SMisPCA=SeqFMisPCA(Sigma,Options);    
            TimeSeqPCA(k,i,l)=toc;
            CorrSeqPCA(k,i,l)=MinFactorsSubspaceDistance(Fo,F_SMisPCA);
            
        end
    end
end

save(['Results/MisPCA_vs_All_vs_dmax_F',num2str(F),'_.mat'])

%% Find SNR point at which <h_est,h> >=Threshold and plot
[h,h2]=PlotFMisPCAComparisonResults(Corr,CorrOr,CorrMisPCA,CorrSeqPCA,dmaxVec,F,SNR,Nrndm,'d_{max}');

figure;
subplot(2,2,1)
imagesc(dmaxVec,SNR,squeeze(mean(TimePCA,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('Plain PCA')
subplot(2,2,2);
imagesc(dmaxVec,SNR,squeeze(mean(TimeOrPCA,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('Plain PCA w/ known delays')
subplot(2,2,3);
imagesc(dmaxVec,SNR,squeeze(mean(TimeMisPCA,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('A-MisPCA')
subplot(2,2,4);
imagesc(dmaxVec,SNR,squeeze(mean(TimeSeqPCA,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('Seq-MisPCA')
suptitle('Ellpased Times (tic-toc)')

figure;
subplot(2,2,1)
imagesc(dmaxVec,SNR,squeeze(mean(Corr,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('Plain PCA')
subplot(2,2,2);
imagesc(dmaxVec,SNR,squeeze(mean(CorrOr,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('Plain PCA w/ known delays')
subplot(2,2,3);
imagesc(dmaxVec,SNR,squeeze(mean(CorrMisPCA,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('A-MisPCA')
subplot(2,2,4);
imagesc(dmaxVec,SNR,squeeze(mean(CorrSeqPCA,2)))
ylabel('SNR')
xlabel('delay magnitude')
title('Seq-MisPCA')
suptitle('Correlation to truth')


saveas(h,['Figures/MisPCA_vs_All_vs_dmax',date,'_F',num2str(F),'_1.fig'],'fig');
saveas(h2,['Figures/MisPCA_vs_All_vs_dmax',date,'_F',num2str(F),'_2.fig'],'fig');

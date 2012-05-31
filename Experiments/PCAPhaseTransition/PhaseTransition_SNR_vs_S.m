%% Sample max eigenval VS Toeplitz bounds, random delays.
addpath(genpath('~/Matlab lib'))
addpath(genpath('../../Utilities'))
clear all;


Nrndm=6*5;

if matlabpool('size') == 0
    matlabpool open 6;
end

Spacing=1;
p=300;
r=round(p/10);
dmax=round(p/10);

NpointsSNR=10;
NpointsN=10;
SNR=round(linspace(-10,30,NpointsSNR));
nVec=round(linspace(100,1e3,NpointsN));
dmaxVec=round(linspace(2,p-1,NpointsN));


%% rectangular
h1=abs([ones(r,1);zeros(p-r,1)]);
h1=h1/norm(h1);

%% rectangular convoluted
h2=abs([ones(r,1);zeros(p-r,1)]);
h2=smoothSig(h2,50,'triang');
h2=h2/norm(h2);

%% sinsusoids
a1=1;a2=1/2;
w1=2*pi*3; w2=2*pi*2;
x=1:p;
h3=(a1*sin(w1*x/p) +a2*sin(w2*x/p))';
h3=h3/norm(h3);

%% Initialization

Corr1=zeros(NpointsSNR,Nrndm,NpointsN);
Corr2=zeros(NpointsSNR,Nrndm,NpointsN);
Corr3=zeros(NpointsSNR,Nrndm,NpointsN);
Corr4=zeros(NpointsSNR,Nrndm,NpointsN);
Corr5=zeros(NpointsSNR,Nrndm,NpointsN);
Corr6=zeros(NpointsSNR,Nrndm,NpointsN);

PT1=zeros(NpointsSNR,Nrndm,NpointsN);
PT2=zeros(NpointsSNR,Nrndm,NpointsN);
PT3=zeros(NpointsSNR,Nrndm,NpointsN);
PT4=zeros(NpointsSNR,Nrndm,NpointsN);
PT5=zeros(NpointsSNR,Nrndm,NpointsN);
PT6=zeros(NpointsSNR,Nrndm,NpointsN);

PT21=zeros(NpointsSNR,Nrndm,NpointsN);
PT22=zeros(NpointsSNR,Nrndm,NpointsN);
PT23=zeros(NpointsSNR,Nrndm,NpointsN);
PT24=zeros(NpointsSNR,Nrndm,NpointsN);
PT25=zeros(NpointsSNR,Nrndm,NpointsN);
PT26=zeros(NpointsSNR,Nrndm,NpointsN);

PTUB1=zeros(NpointsSNR,Nrndm,NpointsN);
PTUB2=zeros(NpointsSNR,Nrndm,NpointsN);
PTUB3=zeros(NpointsSNR,Nrndm,NpointsN);
PTUB4=zeros(NpointsSNR,Nrndm,NpointsN);
PTUB5=zeros(NpointsSNR,Nrndm,NpointsN);
PTUB6=zeros(NpointsSNR,Nrndm,NpointsN);

PTLB1=zeros(NpointsSNR,Nrndm,NpointsN);
PTLB2=zeros(NpointsSNR,Nrndm,NpointsN);
PTLB3=zeros(NpointsSNR,Nrndm,NpointsN);
PTLB4=zeros(NpointsSNR,Nrndm,NpointsN);
PTLB5=zeros(NpointsSNR,Nrndm,NpointsN);
PTLB6=zeros(NpointsSNR,Nrndm,NpointsN);

%% Experiments vs n number of samples)

for l=1:NpointsN
    n=nVec(l);
     fprintf('\n n=%d \n',nVec(l))
    for k=1:NpointsSNR
        fprintf('SNR=%2.2g ',SNR(k))
        % Parallelizable for (deactivated for OCTAVE compatibility) 
 for i=1:Nrndm
            theta=sqrt(10^(SNR(k)/10));
            %% Filtered rectangular
            [S,d,s_d,V]=GenerateMisalignedRank1PlusNoise(theta*h1,n,dmax,Spacing);           
            hest1=eigs(S,1);
            xcorrhhest1=xcorr(V,hest1);
            Corr1(k,i,l)=max(abs(xcorrhhest1));
            [dummy,PT1(k,i,l),PT21(k,i,l),PTUB1(k,i,l),PTLB1(k,i,l)]=PhaseTransitionBounds(h1,d,s_d,1);
            
            %% Filtered rect * triangular
            [S,d,s_d,V]=GenerateMisalignedRank1PlusNoise(theta*h2,n,dmax,Spacing);         
            hest2=eigs(S,1);
            xcorrhhest2=xcorr(V,hest2);
            Corr2(k,i,l)=max(abs(xcorrhhest2));
            [dummy,PT2(k,i,l),PT22(k,i,l),PTUB2(k,i,l),PTLB2(k,i,l)]=PhaseTransitionBounds(h2,d,s_d,1);
                        
            %% Sinusoids
            [S,d,s_d,V]=GenerateMisalignedRank1PlusNoise(theta*h3,n,dmax,Spacing);        
            hest3=eigs(S,1);
            xcorrhhest3=xcorr(V,hest3);
            Corr3(k,i,l)=max(abs(xcorrhhest3));
            [dummy,PT3(k,i,l),PT23(k,i,l),PTUB3(k,i,l),PTLB3(k,i,l)]=PhaseTransitionBounds(h3,d,s_d,1);          
        end
    end
end

%% Experiments vs delay amplitudes
n=p;
for l=1:NpointsN
    dmax=dmaxVec(l);
    fprintf('\n dmax=%d \n',dmax)
    for k=1:NpointsSNR
        fprintf(' SNR=%2.2g ',SNR(k))
        % Parallelizable for (deactivated for OCTAVE compatibility) 
 for i=1:Nrndm
            theta=sqrt(10^(SNR(k)/10));
            %% Filtered rectangular
            [S,d,s_d,V]=GenerateMisalignedRank1PlusNoise(theta*h1,n,dmax,Spacing);           
            hest1=eigs(S,1);
            xcorrhhest1=xcorr(V,hest1);
            Corr4(k,i,l)=max(abs(xcorrhhest1));
            [dummy,PT4(k,i,l),PT24(k,i,l),PTUB4(k,i,l),PTLB4(k,i,l)]=PhaseTransitionBounds(h1,d,s_d,1);
            
            %% Filtered rect * triangular
            [S,d,s_d,V]=GenerateMisalignedRank1PlusNoise(theta*h2,n,dmax,Spacing);         
            hest2=eigs(S,1);
            xcorrhhest2=xcorr(V,hest2);
            Corr5(k,i,l)=max(abs(xcorrhhest2));
            [dummy,PT5(k,i,l),PT25(k,i,l),PTUB5(k,i,l),PTLB5(k,i,l)]=PhaseTransitionBounds(h2,d,s_d,1);
                        
            %% Sinusoids
            [S,d,s_d,V]=GenerateMisalignedRank1PlusNoise(theta*h3,n,dmax,Spacing);        
            hest3=eigs(S,1);
            xcorrhhest3=xcorr(V,hest3);
            Corr6(k,i,l)=max(abs(xcorrhhest3));
            [dummy,PT6(k,i,l),PT26(k,i,l),PTUB6(k,i,l),PTLB6(k,i,l)]=PhaseTransitionBounds(h3,d,s_d,1);          
        end
    end
end


%% Plot results vs n number of samples)
h=figure;
subaxis(3,3,1,'Spacing',0.05,'MR',.1,'ML',0.05); 
plot(h1,'LineWidth',2);
xlabel('Time','FontSize',12)
ylabel('Magnitude','FontSize',12)
title('Rectangular','FontSize',12)

subaxis(3,3,2,'Spacing',0.05,'MR',0.05,'ML',0.05); 
plot(h2,'LineWidth',2);
xlabel('Time','FontSize',12)
title('Rectangular convoluted with triangular','FontSize',12)

subaxis(3,3,3,'Spacing',0.05,'MR',0.05,'ML',0.05); 
plot(h3,'LineWidth',2);
xlabel('Time','FontSize',12)
title('2 sinusoids','FontSize',12)

subaxis(3,3,4,'Spacing',0.05,'MR',.1,'ML',0.05); 
surf(nVec,SNR,squeeze(mean(Corr1.^2,2))); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PT1,2),1))),ones(NpointsN,1),'.-k','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PT21,2),1))),ones(NpointsN,1),'-.r','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PTUB1,2),1))),ones(NpointsN,1),'--^w','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PTLB1,2),1))),ones(NpointsN,1),'--vw','LineWidth',2); hold on;
axis([min(nVec) max(nVec) min(SNR) max(SNR) 0 1])
shading('interp');
view(0,90);
xlabel('n','FontSize',12)
ylabel('SNR','FontSize',12)

subaxis(3,3,5,'Spacing',0.05,'MR',0.05,'ML',0.05); 
surf(nVec,SNR,squeeze(mean(Corr2.^2,2))); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PT2,2),1))),ones(NpointsN,1),'.-k','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PT22,2),1))),ones(NpointsN,1),'.-r','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PTUB2,2),1))),ones(NpointsN,1),'--^w','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PTLB2,2),1))),ones(NpointsN,1),'--vw','LineWidth',2); hold on;
axis([min(nVec) max(nVec) min(SNR) max(SNR) 0 1])
caxis(rangeC)
shading('interp');
view(0,90);
xlabel('n','FontSize',12)

subaxis(3,3,6,'Spacing',0.05,'MR',0.05,'ML',0.05); 
surf(nVec,SNR,squeeze(mean(Corr3.^2,2))); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PT3,2),1))),ones(NpointsN,1),'.-k','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PT23,2),1))),ones(NpointsN,1),'.-r','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PTUB3,2),1))),ones(NpointsN,1),'--^w','LineWidth',2); hold on;
plot3(nVec,10*log10(squeeze(mean(mean(PTLB3,2),1))),ones(NpointsN,1),'--vw','LineWidth',2); hold on;
axis([min(nVec) max(nVec) min(SNR) max(SNR) 0 1])
caxis(rangeC)
shading('interp');
view(0,90);
xlabel('n','FontSize',12)

%% Plot results vs dmax
subaxis(3,3,7,'Spacing',0.05,'MR',.1,'ML',0.05,'MarginBottom',.1); 
surf(dmaxVec,SNR,squeeze(mean(Corr4.^2,2))); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PT4,2),1))),ones(NpointsN,1),'.-k','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PT24,2),1))),ones(NpointsN,1),'-.r','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PTUB4,2),1))),ones(NpointsN,1),'--^w','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PTLB4,2),1))),ones(NpointsN,1),'--vw','LineWidth',2); hold on;
axis([min(dmaxVec) max(dmaxVec) min(SNR) max(SNR) 0 1])
caxis(rangeC)
shading('interp');
view(0,90);
ylabel('SNR','FontSize',12)
xlabel('d_{max}','FontSize',12)

subaxis(3,3,8,'Spacing',0.05,'MR',0.05,'ML',0.05,'MarginBottom',.1); 
surf(dmaxVec,SNR,squeeze(mean(Corr5.^2,2))); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PT5,2),1))),ones(NpointsN,1),'.-k','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PT25,2),1))),ones(NpointsN,1),'.-r','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PTUB5,2),1))),ones(NpointsN,1),'--^w','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PTLB5,2),1))),ones(NpointsN,1),'--vw','LineWidth',2); hold on;
axis([min(dmaxVec) max(dmaxVec) min(SNR) max(SNR) 0 1])
caxis(rangeC)
shading('interp');
view(0,90);
xlabel('d_{max}','FontSize',12)

subaxis(3,3,9,'Spacing',0.05,'MR',0.05,'ML',0.05,'MarginBottom',.1); 
surf(dmaxVec,SNR,squeeze(mean(Corr6.^2,2))); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PT6,2),1))),ones(NpointsN,1),'.-k','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PT26,2),1))),ones(NpointsN,1),'.-r','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PTUB6,2),1))),ones(NpointsN,1),'--^w','LineWidth',2); hold on;
plot3(dmaxVec,10*log10(squeeze(mean(mean(PTLB6,2),1))),ones(NpointsN,1),'--vw','LineWidth',2); hold on;
axis([min(dmaxVec) max(dmaxVec) min(SNR) max(SNR) 0 1])
caxis(rangeC)
shading('interp');
view(0,90);
xlabel('d_{max}','FontSize',12)
colorbar('SouthOutside')


saveas(h,['Figures/PT_PCA_',date,'.fig'],'fig');
saveas(h,['Figures/PT_PCA_',date,'.eps'],'eps');

save(['Results/PT_PCA_',date,'.mat'])

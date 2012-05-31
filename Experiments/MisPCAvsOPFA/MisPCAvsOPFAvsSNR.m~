close all
clear all

N_CPU=4;

if matlabpool('size') == 0
    matlabpool('open', N_CPU) ;
end

addpath(genpath('../../Utilities'));
addpath(genpath('../../../OPFA_package/Utilities'));

p=50;
f_o=2;
f_no=0;
n=50;
n_f=n;
S=10;
Nrndm=N_CPU*5;
VarDelays=5;
dmax=round(sqrt(12*VarDelays+1));
FigOnOff=0;
Npoints=12;
SNRVec=linspace(-30,30,Npoints);

Err_OPFA=zeros(Npoints,Nrndm);
Err_FMisPCA=zeros(Npoints,Nrndm);
MErr_OPFA=zeros(Npoints,Nrndm);
MErr_FMisPCA=zeros(Npoints,Nrndm);
MErr_PF=zeros(Npoints,Nrndm);
MErr_PF2=zeros(Npoints,Nrndm);
RErr_OPFA=zeros(Npoints,Nrndm);
RErr_DL=zeros(Npoints,Nrndm);
RErr_PF=zeros(Npoints,Nrndm);
RErr_PF2=zeros(Npoints,Nrndm);
SNRVecEst=zeros(Npoints,Nrndm);
Rank_FMisPCA=zeros(Npoints,Nrndm);

for i=1:Npoints
    SNR=SNRVec(i);
    fprintf('*************** SNR=%d (%d/%d)\n',SNR,i,Npoints)
    Aux1=zeros(1,Nrndm);
    Aux2=zeros(1,Nrndm);
    Aux3=zeros(1,Nrndm);
    Aux4=zeros(1,Nrndm);
    Aux11=zeros(1,Nrndm);
    Aux12=zeros(1,Nrndm);
    Aux13=zeros(1,Nrndm);
    Aux14=zeros(1,Nrndm);
    ErrAux1=zeros(1,Nrndm);
    ErrAux2=zeros(1,Nrndm);
    ErrAux3=zeros(1,Nrndm);
    ErrAux4=zeros(1,Nrndm);

    SNRAux=zeros(1,Nrndm);
    
    % Parallelizable for (deactivated for OCTAVE compatibility) 
 for j=1:Nrndm

        fprintf('******************** MC %d/%d \n',j,Nrndm)
        [A_o,U,F,X,D,dini,SNRAux(j)]=GenerateFMisPCADataForOPFA(p,n,n_f,S,f_o,SNR,VarDelays,0);        
        [Options]=OPFA_Options(X,f_o,0,dmax);
        [Options.LambdaMax]=EstimateLambdaMax(X,Options.Masks,Options);
        Options.Initialization=1;
        Options.Nrndm=4;
        
        disp('OPFA')

        Options.BetaVec=logspace(-6,-3,1);
        Options.LambdaVec=logspace(-8,log10(Options.LambdaMax),5);
        [OutputOPFA]=OPFA(X,Options);
        
        disp('MisPCA')
       
        FMisPCAOptions=FMisPCA_options([],n,p,1);
        MisPCAResults=FMisPCAwCV(X,FMisPCAOptions);

        [dummy1,dummy2,Aux1(j),Aux11(j)]=OPFA_Performance(X,D,OutputOPFA{1},F,dini,FigOnOff,'OPFA');
        [dummy1,dummy2,Aux2(j),Aux12(j)]=FMisPCA_Performance(X,D,Options.Masks,F,MisPCAResults.F,dini,MisPCAResults.d,MisPCAResults.Lambdas,FigOnOff);
        RankAux(j)=MisPCAResults.f;
    end
    for j=1:Nrndm
        Err_OPFA(i,j)=Aux1(j);
        Err_FMisPCA(i,j)=Aux2(j);

        MErr_OPFA(i,j)=Aux11(j);
        MErr_FMisPCA(i,j)=Aux12(j);
        
        Rank_FMisPCA(i,j)=RankAux(j);
        
        SNRVecEst(i,j)=SNRAux(j);
    end
    %% Stock stats
    Stats{1}=Err_OPFA;
    Stats{2}=Err_FMisPCA;
    
    MStats{1}=MErr_OPFA;
    MStats{2}=MErr_FMisPCA;
    
    xstr=mean(SNRVecEst,2);
    legendStr={'OPFA','SFA','OPFA-C'};
    Results.xstr=xstr;
    Results.Stats=Stats;
    Results.MStats=MStats;
    Results.Rank_FMisPCA=Rank_FMisPCA;
    Results.legendStr=legendStr;
    
    save(['Results/MisPCAvsOPFA',strrep(date,'-','_'),'.mat'],'Results')
    
%     h=figure;
%     subplot(2,1,1)
%     plot(xstr,mean(Stats{1},2),'-or'); hold on;
%     plot(xstr,mean(Stats{2},2),'-xb'); hold on;
%     ylabel('Distance between subspaces','FontSize',12)
%     leg=legend('OPFA','MisPCA');
%     set(leg,'FontSize',12);
%     grid on;
%     subplot(2,1,2)
%     semilogy(xstr,mean(MStats{1},2),'-or'); hold on;
%     semilogy(xstr,mean(MStats{2},2),'-xb'); hold on;
%     xlabel('SNR','FontSize',12)
%     ylabel('MSE','FontSize',12)
%     grid on;
%     if i< Npoints
%         saveas(h,strcat('Figures/OPFAvsMisPCA',strrep(date,'-','_'),'_i=',num2str(i),'SNR.fig'))
%     end
end

h=figure;
subplot(2,1,1)
plot(xstr,mean(Stats{1},2),'-or'); hold on;
plot(xstr,mean(Stats{2},2),'-xb'); hold on;
ylabel('Distance between subspaces','FontSize',12)
leg=legend('OPFA','MisPCA');
set(leg,'FontSize',12);
grid on;
subplot(2,1,2)
semilogy(xstr,mean(MStats{1},2),'-or'); hold on;
semilogy(xstr,mean(MStats{2},2),'-xb'); hold on;
xlabel('SNR','FontSize',12)
ylabel('MSE','FontSize',12)
grid on;
saveas(h,strcat('Figures/MisPCAvsOPFA',strrep(date,'-','_'),'_SNR.fig'))
suptitle('MisPCA model')

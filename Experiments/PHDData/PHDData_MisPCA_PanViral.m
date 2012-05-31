close all
clear all

addpath(genpath('../../Utilities'));

N_CPU=6;
if matlabpool('size') == 0
   matlabpool('open', N_CPU) ;
end

%% Choose Data
p=300;
Padding=4;
Options.Padding=Padding;
FixedDelaysOnOff=1; % Fix all factor's delays to be the same or not.
InterpOnOff=1;

%% Load Data and Gene List
Virus='C';
OldOrNew='New';
SxOrAsx=1;

[sample,Symptoms,TimeVecOr1]=LoadData('./Data/',Virus,OldOrNew,0);
TimeVecOr1(1)=-12;
AllGenes=sample.rownames(1:end-62);

[Options,OnsetTimes1,Labels]=LoadSymptomsTable(SxOrAsx,Virus,OldOrNew,Symptoms,Options);

Options.SubsMeanOnOff=1;
[X1,Masks1,AllTimes,Samples,y1,RowIds,AllGenes,GenesAffy]=ExtractSubjectsRNAProfile(sample,Options.SxSubjectIDs,...
                                                    Virus,AllGenes,0,0,Options.SubsMeanOnOff,0);                                  
Virus='W';
OldOrNew='New';
[sample,Symptoms2,TimeVecOr2]=LoadData('./Data/',Virus,OldOrNew,0);
TimeVecOr2(1)=-12;
[Options2,OnsetTimes2,Labels2]=LoadSymptomsTable(SxOrAsx,Virus,OldOrNew,Symptoms2,Options);
[X2,Masks2,~,~,y2]=ExtractSubjectsRNAProfile(sample,Options2.SxSubjectIDs,...
                                                    Virus,AllGenes,0,0,Options.SubsMeanOnOff,0);
                                                
                                              
Virus='Z';
OldOrNew='Old';
[sample,Symptoms3,TimeVecOr3]=LoadData('./Data/',Virus,OldOrNew,0);
[Options3,OnsetTimes3,Labels3]=LoadSymptomsTable(SxOrAsx,Virus,OldOrNew,Symptoms3,Options);

[X3,Masks3,~,~,y3]=ExtractSubjectsRNAProfile(sample,Options3.SxSubjectIDs,...
                                                    Virus,AllGenes,0,0,Options.SubsMeanOnOff,0);
                                                
Xc={X1,X2,X3};
Masksc={Masks1,Masks2,Masks3};
yc={y1,y2,y3};
TimeVecOrc={TimeVecOr1,TimeVecOr2,TimeVecOr3};
OnsetTimes=[OnsetTimes1;OnsetTimes2;OnsetTimes3];

[X,Masks,y,TimeVecOr]=combine_PHD_datasets(Xc,Masksc,TimeVecOrc,yc,Options.Padding);
Options.n=size(X{1},1);
Options.S=length(X);
Options.p=p;

Times=[1:4';5:8';9:12';13:16']; % Time points to compare in the ANOVA procedure

[Genes,GenesAffy,RowIds,X,Masks,SortedPv]=SelectTimeCourseGenesANOVA(AllGenes,GenesAffy,RowIds,Times,X,Masks,Options);

if InterpOnOff
    [X,Masks]=InterpolateMissing(X,Masks);
end

%% Initialization
[Options]=PHD_options(Options,p,Options.n,N_CPU,FixedDelaysOnOff);
Options.kvec=1:4;

[Results]=MisPCAwCV(X,Masks,Options);

PlotMisPCAwCVResults(Results);
    
[Err,h]=OriginalVSReconstructed(Results.d,Results.F,X,X,Masks,Options.Padding,...
    [1:100],[1:Options.S],1,[1:Options.S],'',Genes,OnsetTimes,Results.Lambdas,TimeVecOr);

[PCAResults]=PCAwCV(X,Masks,Options)

%% Save results 

Results.Options=Options;
Results.Data.Genes=Genes;
Results.Data.GenesAffy=GenesAffy;
Results.Data.GenesPv=SortedPv;
Results.Data.X=X;
Results.Data.Masks=Masks;
Results.Data.y=y;
Results.Data.TimeVecOr=TimeVecOr;
Results.Data.Labels=Labels;
Results.Data.OnsetTimes=OnsetTimes;

%Results.PCA=PCAResults;

if SxOrAsx==1
  TypeFile='Sx';
else
  TypeFile='Asx';
end

save(strcat('Results/MisPCA_PanViral_',TypeFile,'_',strrep(date,'-','_'),'_p=',num2str(Options.p),'_',Virus,'_',OldOrNew,'.mat'), 'Results');
%% Plot  results

[Err,h]=OriginalVSReconstructed(PCAResults.d,PCAResults.F,X,X,Masks,Options.Padding,...
    [1:100],[1:Options.S],1,[1:Options.S],'',Genes,OnsetTimes,PCAResults.Lambdas,TimeVecOr);


%% Cluster Analysis
%% Variable Selection and Cluster Analysis
newp=500;
[A_misPCA,A_PCA,Xvs,Genesvs,GenesAffyvs,SortedPvvs]=ConcatLoadings(Results,PCAResults,Options,X,Genes,GenesAffy,SortedPv,newp);


linkage='single';
distance='seuclidean';
Nclust=8

TMisPCA = clusterdata(A_misPCA,'maxclust',Nclust,'distance',distance,'linkage',linkage); 
TPCA = clusterdata(A_PCA,'maxclust',Nclust,'distance',distance,'linkage',linkage); 

% [t,in]=sort(TMisPCA);
% figure;imagesc(A_misPCA(in,:));
% 
GenesForClustering=1:newp;
ClusterLoadingsOrResponses=2;
Showclus=[1:Nclust];

[Signature,ClusterMembersMisPCA,XrMisPCA,AlignedDataMisPCA,T1]=TemporalSignatures(TMisPCA,Results.F,Xvs,Masks,Results.d,Results.Lambdas,...
    GenesForClustering,Genesvs,GenesAffyvs,SortedPvvs,Options,Results.Data.TimeVecOr,ClusterLoadingsOrResponses,FixedDelaysOnOff,Showclus);

[SignaturePCA,ClusterMembersPCA,XrPCA,AlignedData,T2]=TemporalSignatures(TPCA,PCAResults.F,Xvs,Masks,PCAResults.d,PCAResults.Lambdas,...
    GenesForClustering,Genesvs,GenesAffyvs,SortedPvvs,Options,Results.Data.TimeVecOr,ClusterLoadingsOrResponses,FixedDelaysOnOff,Showclus);

[Corr]=FindSignatureCorrespondance(Signature,SignaturePCA)  

OrderClustersMisPCA=[1,3,4,2,6,7]; %[6,5,4,3,1,2];           
OrderClustersPCA=[1,2,3,8,6,5]; %1,2 

Plot_MisPCAClusters(TMisPCA,TPCA,Signature,SignaturePCA,Results.A,PCAResults.A,TimeVecOr,OrderClustersMisPCA,OrderClustersPCA,Options)            


SelectedSubjects=[1,4,10,12,14,16,18,24]
LabelsSubjects={'HRV','HRV','HRV','H1N1','H1N1','H3N2','H3N2','H3N2'};
figure;
for i=1:length(SelectedSubjects)
    ni=length(find(Masks{SelectedSubjects(i)}(:,1)==1));
    TimeIndx=Options.Padding+1:(ni-Options.Padding);
    Indxtime=(Masks{SelectedSubjects(i)}(:,1)==1);
    subplot(length(SelectedSubjects),1,i)
    plot(TimeVecOr(Indxtime),Xvs{SelectedSubjects(i)});
    axis tight;
    title(['Subject ID:',num2str(SelectedSubjects(i)),'(',LabelsSubjects{i},')'])
    ylabel('Gene Expression','FontSize',12)
    if i==length(SelectedSubjects)
        xlabel('Time (Hours)','FontSize',12)
    end
end
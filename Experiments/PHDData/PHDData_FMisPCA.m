close all;
clear all;
addpath(genpath('../../Utilities'));

N_CPU=6;
if matlabpool('size') == 0
   matlabpool('open', N_CPU) ;
end

%% Choose Data
SxOrAsx=1; % '0' Asx, '1' Sx, '2' Other
p=1000;
Padding=4;
Options.Padding=Padding;
FixedDelaysOnOff=1; % Fix all factor's delays to be the same or not.
InterpOnOff=0;

%% Load Data and Gene List
Virus='Z';
OldOrNew='Old';
SxOrAsx=1;
[sample,Symptoms,TimeVecOr]=LoadData('./Data/',Virus,OldOrNew,Options.Padding);

AllGenes=sample.rownames(1:end-62);

[Options,OnsetTimes,Labels]=LoadSymptomsTable(SxOrAsx,Virus,OldOrNew,Symptoms,Options);

Options.SubsMeanOnOff=1;
[X,Masks,AllTimes,Samples,y,RowIds,AllGenes,GenesAffy]=ExtractSubjectsRNAProfile(sample,Options.SxSubjectIDs,...
                                                    Virus,AllGenes,0,Options.Padding,Options.SubsMeanOnOff,InterpOnOff);

Times=[1:4';5:8';9:12';13:16']; % Time points to compare in the ANOVA procedure

Options.p=p;
Options.n=size(Masks{1},1);
[Genes,GenesAffy,RowIds,X,Masks,SortedPv]=SelectTimeCourseGenesANOVA(AllGenes,GenesAffy,RowIds,Times,X,Masks,Options);

Options.n=size(Masks{1},1);


%% Initialization
[Options]=PHD_options(Options,p,Options.n,N_CPU,FixedDelaysOnOff);
Options.dmax=10;

[Results]=FMisPCAwCV(X,Options)

PlotMisPCAwCVResults(Results);

[Err,h]=OriginalVSReconstructed(Results.d,Results.F,X,X,Masks,Options.Padding,...
    [1:100],[1:Options.S],1,Options.SxSubjectIDs,'',Genes,OnsetTimes,Results.Lambdas,TimeVecOr);


%break
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

save(strcat('Results/OPPCA_only',TypeFile,'_',strrep(date,'-','_'),'_p=',num2str(Options.p),'_',Virus,'_',OldOrNew,'.mat'), 'Results','PCAResults');
%% Plot  results

[Err,h]=OriginalVSReconstructed(PCAResults.d,PCAResults.F,X,X,Masks,Options.Padding,...
    [1:100],[1:Options.S],1,Options.SxSubjectIDs,'',Genes,OnsetTimes,PCAResults.Lambdas,TimeVecOr);

%% Variable Selection and Cluster Analysis
newp=1000;
AverageOnOff=1;
[A_misPCA,A_PCA,Xvs,Genesvs,GenesAffyvs,SortedPvvs]=ConcatLoadings(Results,PCAResults,Options,X,Genes,GenesAffy,SortedPv,newp,AverageOnOff);

linkage='complete';
distance='seuclidean';
Nclust=6

TMisPCA = clusterdata(A_misPCA,'maxclust',Nclust,'distance',distance,'linkage',linkage); 
TPCA = clusterdata(A_PCA,'maxclust',Nclust,'distance',distance,'linkage',linkage); 
Nclust=min(length(unique(TMisPCA)),length(unique(TPCA)));

[t,in]=sort(TMisPCA);
figure;
imagesc([repmat(t(in),1,3),A_misPCA(in,:)]);

GenesForClustering=1:newp;
ClusterLoadingsOrResponses=1;
Showclus=[1:Nclust];

[Signature,ClusterMembersMisPCA,XrMisPCA,AlignedDataMisPCA,T1]=TemporalSignatures(TMisPCA,Results.F,Xvs,Masks,Results.d,Results.Lambdas,...
    GenesForClustering,Genesvs,GenesAffyvs,SortedPvvs,Options,Results.Data.TimeVecOr,ClusterLoadingsOrResponses,FixedDelaysOnOff,Showclus);

[SignaturePCA,ClusterMembersPCA,XrPCA,AlignedData,T2]=TemporalSignatures(TPCA,PCAResults.F,Xvs,Masks,PCAResults.d,PCAResults.Lambdas,...
    GenesForClustering,Genesvs,GenesAffyvs,SortedPvvs,Options,Results.Data.TimeVecOr,ClusterLoadingsOrResponses,FixedDelaysOnOff,Showclus);

[Corr]=FindSignatureCorrespondance(Signature,SignaturePCA);

OrderClustersMisPCA=[1,2,6,3,4,5]; %[6,5,4,3,1,2];           
OrderClustersPCA=[Corr(OrderClustersMisPCA)]; %1,2 

h=Plot_MisPCAClusters(TMisPCA,TPCA,Signature,SignaturePCA,A_misPCA,A_PCA,TimeVecOr,OrderClustersMisPCA,OrderClustersPCA,ClusterMembersMisPCA,ClusterMembersPCA,Options)

saveas(h,strcat('figures/OPPCA_only',TypeFile,'_',strrep(date,'-','_'),'_p=',num2str(Options.p),'_',Virus,'_',OldOrNew,'.fig')')
            
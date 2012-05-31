function [Genes]=TimeCourseGenesList(DataDir,sample,p)

if p<length(sample.rownames)
    Genes={'HOXA6','CD177','CD72','OAS1','OASL','ADA','CCR1','LY6E','CD1C','ORM1','RPL3','CX3CR1','ALDH5A1','LY6E', 'CCRL2', 'CCR5','CCL5'};
    % GenesAsx=load(strcat(DataDir,'x100GenesAsxInf.mat'));
    % [GenesAsx,GeneListaffynew]=Convertaffy2ccdf(GenesAsx.x100GenesAsxInf,DataDir);
    % GenesAsx=GenesAsx([1:50 52:100]);
    GenesSx=load(strcat(DataDir,'x568DEGenes.mat'));
    % GenesSx=load(strcat(DataDir,'x900DEGenes.mat'));
    % GenesSx.DEGenes.Genes=GenesSx.DEGenes.Genes(3:length(GenesSx.DEGenes.Genes));
    % Genes=union(Genes,GenesAsx);
    Genes=union(Genes,GenesSx.DEGenes.Genes(1:min(p,length(GenesSx.DEGenes.Genes))-length(Genes)));
    length(Genes);
    if length(Genes)<p
        complGenes=sample.rownames(63:end);
        LeftGenes=setdiff(complGenes,Genes);
        LeftGenes=LeftGenes(20:length(LeftGenes));
        Genes=union(Genes,LeftGenes(1:(p-length(Genes))));
    end
else
    Genes=sample.rownames;
end

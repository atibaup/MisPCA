function [sample,Symptoms,TimeVecOr]=LoadData(DataDir,Virus,OldOrNew,Padding)

if strcmp(Virus,'Z')
    load(strcat(DataDir,'datastructures_RNA_H3N2_interp'));
elseif strcmp(Virus,'W')
    if strcmp(OldOrNew,'Old')
        load(strcat(DataDir,'datastructures_RNA_H1N1_interp'));
    else
        load(strcat(DataDir,'datastructures_H1N1'));
    end
elseif strcmp(Virus,'C')
    load(strcat(DataDir,'datastructures_HRV'));
end
load(strcat(DataDir,'Symptoms'));
sample.blood.RNA=sample;
TimeVecOr=[(-Padding:-1)+min(sample.excelmat_clock(1,:)),sample.excelmat_clock(1,:),max(sample.excelmat_clock(1,:))+(1:Padding)];

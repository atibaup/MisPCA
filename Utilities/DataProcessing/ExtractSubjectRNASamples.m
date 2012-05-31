function [Data]=ExtractSubjectRNASamples(sample,subjectID,Virus)

[Data.Rows]=FindSubjectRows(num2str(subjectID),Virus,sample.blood.RNA.subject);
Data.Samples=intersect(find(sample.blood.RNA.colnames(:,1)==Virus),Data.Rows);
[order reindx]=sort(sample.blood.RNA.time(Data.Samples));
Data.subjectID=subjectID;
Data.Samples=Data.Samples(reindx);

Data.Subjects=sample.blood.RNA.subject(Data.Samples);
Data.SampleTimes=sample.blood.RNA.time(Data.Samples);
Data.colnames=sample.blood.RNA.colnames(Data.Samples,:);
Data.Label=sample.blood.RNA.pheno(Data.Samples,:);
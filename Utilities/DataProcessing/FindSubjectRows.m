function [Rows]=FindSubjectRows(Subject,Virus,SubjectRowsCell)
if length(Subject)==1
    Subject=strcat('0',Subject);
end
VirSub=strcat(Virus,Subject);
k=1;
Rows=[];
for i=1:length(SubjectRowsCell)
    if iscell(SubjectRowsCell)
        if strcmp(SubjectRowsCell{i},VirSub)
            Rows(k)=i;
            k=k+1;
        end
    else
        if strcmp(SubjectRowsCell(i),VirSub)
            Rows(k)=i;
            k=k+1;
        end
    end
end
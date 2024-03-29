function samplmtrx=smpl_select(rowsel,colsel,excelmat,imtrx,alltimes,diagnostic);
%function samplmtrx=smpl_select(excelvec,imtrx);
%Selector function for translating sampling times into chip samples
%Inputs
%rowsel - imtrx rows selected
%colsel - imtrx columns selected
%excelvec - sampling times 
%imtrx - subj be time chip index matrix generated by 
%alltimes - vector of times (length equals number cols of imtrx)
%diagnostic - 1 if missing chips are to be flagged, 0 otherwise
%Outputs
%samplmtrx - the vector of corresponding chips (zero entry means missingchip)
%Example: to generate 0.2T chips for all sx subjects: samplmtrx=smpl_select(find(Subjpheno==0),4,excelmat,imtrx,alltimes)
%-----------------------------------
samplmtrx=zeros(length(rowsel),length(colsel));
for i=1:length(rowsel)
    for j=1:length(colsel)
        if isempty(find(alltimes==excelmat(rowsel(i),colsel(j))))
if diagnostic==1
    ['Bad selected time i=',num2str(i),', j=',num2str(j),' entry set to zero']
end
        else
%            find(alltimes==excelmat(rowsel(i),colsel(j)))
          samplmtrx(i,j)=imtrx(rowsel(i),find(alltimes==excelmat(rowsel(i),colsel(j))));
        end
    end
end
end
function pv=pv_behrens_fisher2(data,timesel,comparison_method,imtrx,excelmat,Subjpheno,alltimes,diagnostic)
%function [pv,pvinds]=pv_behrens_fisher(data,timesel,comparison_method,imtrx,excelmat,Subjpheno,alltimes)
% Uses the critical region of the unpaired Behren's Fisher t-test for computing pvalues 
% on the difference between means of two groups of columns A and B of "data", specified by the indicator 
% "Subjpheno,"  over a time course specified by "timessel." 
% A. Hero sept 2008.
%-----------------
% Inputs
% data: matrix of data (rows are genes, columns are samples)
% timesel: times to aggregate the chips, out of 1,2,3,4,5,6 
% comparison_method = 
%                   all, all absolute times (acording to timesel) are aggregated together
%                   allT, all postchallenge symptom times (0.1,0.2,0.8, T)
%                   maxT, max_t(abs(mean(A_t-S_t)))>0 significant (instantaneous dominance over timesel)
%                   minT, min_t(abs(mean(A_t-S_t)))>0 significant (persistent dominance over timesel)
%                   DeltaTphenomax, same as maxT except A_t-S_t replaced by DeltaA_t-DeltaS_t=(A_t-A_{t-1})-(S_{t}-S_{t-1})
%                   DeltaTphenomin, same as maxT except A_t-S_t replaced by DeltaA_t-DeltaS_t
%                   DeltaTmax, max_t(abs(mean(DeltaA))) OR max_t(abs(mean(DeltaS)))>0 significant 
%                   DeltaTmin, min(min_t(abs(mean(DeltaA))) OR min_t(abs(mean(DeltaS)))>0 significant 
% imtrx: the matrix mapping subject vs time to chip (read_IA_data.m)
% excelmat: the list of BL, PC, 0.1T, etc
% diagnotic: 1 if missing chips for given times is to reported
% Subjpheno: the phenotype of each subject (length is row dimension of imtrx), this is the indicator vector to separate classes  
% alltimes: all times that samples were collects (length is column dimension of imtrx)
% Outputs
% pv: vector of p values
%--------------------------
% Key for timesel
%1    2     3     4     5      6
%-12  0     0.1T  0.2T  0.8T   T
%-------------------------
q=0.2; % necessary for ttest2.m call
[Ngenes,nc]=size(data);
ginds=1:Ngenes; 


switch comparison_method
    case 'all'
        

%Collect all samples from each pheno
%[nr,nc]=size(imtrx(find(Subjpheno==0),timesel));samplmtrx_S=reshape(imtrx(find(Subjpheno==0),timesel),nc*nr,1);
%[nr,nc]=size(imtrx(find(Subjpheno==1),timesel));samplmtrx_A=reshape(imtrx(find(Subjpheno==1),timesel),nc*nr,1);


samplmtrx_oneandtwoT_A=smpl_select(find(Subjpheno==1),timesel,excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrx_oneandtwoT_S=smpl_select(find(Subjpheno==0),timesel,excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
[r,c]=size(samplmtrx_oneandtwoT_A);
samplmtrx_A=reshape(samplmtrx_oneandtwoT_A,r*c,1);
[r,c]=size(samplmtrx_oneandtwoT_S);
samplmtrx_S=reshape(samplmtrx_oneandtwoT_S,r*c,1);
samplmtrx_A(samplmtrx_A==0)=[]; %eliminate missing chips
samplmtrx_S(samplmtrx_S==0)=[];
%colnames(samplmtrx_A,:);
%colnames(samplmtrx_S,:); 
Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);
['[Nsamp_A,Nsamp_S]=',num2str([Nsamp_A,Nsamp_S])]
    
[h,pv]=ttest2(data(ginds,samplmtrx_A)',data(ginds,samplmtrx_S)',q,'both','unequal');

    
    case 'allT'
        

samplmtrx_oneandtwoT_A=smpl_select(find(Subjpheno==1),[3,4,5,6],excelmat,imtrx,alltimes,diagnostic)%Identify samples from 0.8T and T
samplmtrx_oneandtwoT_S=smpl_select(find(Subjpheno==0),[3,4,5,6],excelmat,imtrx,alltimes,diagnostic)%Identify samples from 0.8T and T
[r,c]=size(samplmtrx_oneandtwoT_A);
samplmtrx_A=reshape(samplmtrx_oneandtwoT_A,r*c,1);
[r,c]=size(samplmtrx_oneandtwoT_S);
samplmtrx_S=reshape(samplmtrx_oneandtwoT_S,r*c,1);
samplmtrx_A(samplmtrx_A==0)=[]; %eliminate missing chips
samplmtrx_S(samplmtrx_S==0)=[];

Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);
['[Nsamp_A,Nsamp_S]=',num2str([Nsamp_A,Nsamp_S])]


[h,pv]=ttest2(data(ginds,samplmtrx_A)',data(ginds,samplmtrx_S)',q,'both','unequal');


    case 'maxT'


pvmin=ones(1,Ngenes);
for T=timesel                 
samplmtrxT_A=smpl_select(find(Subjpheno==1),T,excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),T,excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
[r,c]=size(samplmtrxT_A);
samplmtrx_A=reshape(samplmtrxT_A,r*c,1);
[r,c]=size(samplmtrxT_S);
samplmtrx_S=reshape(samplmtrxT_S,r*c,1);
samplmtrx_A(samplmtrx_A==0)=[]; %eliminate missing chips
samplmtrx_S(samplmtrx_S==0)=[];

Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);
['[Nsamp_A,Nsamp_S]=',num2str([Nsamp_A,Nsamp_S])]

[h,pv]=ttest2(data(ginds,samplmtrx_A)',data(ginds,samplmtrx_S)',q,'both','unequal');

pvmin=min(pvmin,pv);
T
end

pv=1-(1-pvmin).^length(timesel);  


    case 'minT'
       


pvmax=zeros(1,Ngenes);
for T=timesel                 
samplmtrxT_A=smpl_select(find(Subjpheno==1),T,excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),T,excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
[r,c]=size(samplmtrxT_A);
samplmtrx_A=reshape(samplmtrxT_A,r*c,1);
[r,c]=size(samplmtrxT_S);
samplmtrx_S=reshape(samplmtrxT_S,r*c,1);
samplmtrx_A(samplmtrx_A==0)=[]; %eliminate missing chips
samplmtrx_S(samplmtrx_S==0)=[];

%samplmtrx_A
%samplmtrx_S %

Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);
['[Nsamp_A,Nsamp_S]=',num2str([Nsamp_A,Nsamp_S])]

[h,pv]=ttest2(data(ginds,samplmtrx_A)',data(ginds,samplmtrx_S)',q,'both','unequal');

pvmax=max(pvmax,pv);
T
end

pv=pvmax.^length(timesel);  

    case 'DeltaTmin'
       


pvmaxA=zeros(1,Ngenes);pvmaxS=zeros(1,Ngenes);
for T=1:length(timesel)-1 %minumum must be greater than 1                 
samplmtrxT_A=smpl_select(find(Subjpheno==1),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
[r,c]=size(samplmtrxT_A);
samplmtrxm_A=reshape(samplmtrxT_A,r*c,1);
[r,c]=size(samplmtrxT_S);
samplmtrxm_S=reshape(samplmtrxT_S,r*c,1);
samplmtrxm_A(samplmtrxm_A==0)=[]; %eliminate missing chips
samplmtrxm_S(samplmtrxm_S==0)=[];

Nsampm_A=length(samplmtrxm_A);
Nsampm_S=length(samplmtrxm_S);

samplmtrxT_A=smpl_select(find(Subjpheno==1),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
[r,c]=size(samplmtrxT_A);
samplmtrx_A=reshape(samplmtrxT_A,r*c,1);
[r,c]=size(samplmtrxT_S);
samplmtrx_S=reshape(samplmtrxT_S,r*c,1);
samplmtrx_A(samplmtrx_A==0)=[]; %eliminate missing chips
samplmtrx_S(samplmtrx_S==0)=[];

Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);



[hA,pvA]=ttest2(data(ginds,samplmtrx_A)',data(ginds,samplmtrxm_A)',q,'both','unequal');
[hS,pvS]=ttest2(data(ginds,samplmtrx_S)',data(ginds,samplmtrxm_S)',q,'both','unequal');
pvmaxA=max(pvmaxA,pvA);
pvmaxS=max(pvmaxS,pvS);
end

%pv=1-(1-(pvmaxA).^(length(timesel)-1)).*(1-(pvmaxS).^(length(timesel)-1));  

pv=min((pvmaxA).^(length(timesel)-1),(pvmaxS).^(length(timesel)-1));  
%pvmaxA(1:10),pvmaxS(1:10),pv(1:10),(length(timesel)-1)

    case 'DeltaTmax'
 

pvminA=ones(1,Ngenes);pvminS=ones(1,Ngenes);
for T=1:length(timesel)-1 %minumum must be greater than 1  
samplmtrxT_A=smpl_select(find(Subjpheno==1),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
[r,c]=size(samplmtrxT_A);
samplmtrxm_A=reshape(samplmtrxT_A,r*c,1);
[r,c]=size(samplmtrxT_S);
samplmtrxm_S=reshape(samplmtrxT_S,r*c,1);
samplmtrxm_A(samplmtrxm_A==0)=[]; %eliminate missing chips
samplmtrxm_S(samplmtrxm_S==0)=[];

Nsampm_A=length(samplmtrxm_A);
Nsampm_S=length(samplmtrxm_S);

samplmtrxT_A=smpl_select(find(Subjpheno==1),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
[r,c]=size(samplmtrxT_A);
samplmtrx_A=reshape(samplmtrxT_A,r*c,1);
[r,c]=size(samplmtrxT_S);
samplmtrx_S=reshape(samplmtrxT_S,r*c,1);
samplmtrx_A(samplmtrx_A==0)=[]; %eliminate missing chips
samplmtrx_S(samplmtrx_S==0)=[];

Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);


[hA,pvA]=ttest2(data(ginds,samplmtrx_A)',data(ginds,samplmtrxm_A)',q,'both','unequal');
[hS,pvS]=ttest2(data(ginds,samplmtrx_S)',data(ginds,samplmtrxm_S)',q,'both','unequal');
pvminA=min(pvminA,pvA);
pvminS=min(pvminS,pvS);
end

%pv=1-(1-(pvmaxA).^(length(timesel)-1)).*(1-(pvmaxS).^(length(timesel)-1));  

pv=min(1-(1-pvminA).^(length(timesel)-1),1-(1-pvminS).^(length(timesel)-1));
%pvminA(1:10),pvminS(1:10),pv(1:10),(length(timesel)-1)

    case 'DeltaTphenomin'
      


pvmax=zeros(1,Ngenes);
for T=1:length(timesel)-1 %minumum must be greater than 1                 
samplmtrxTm_A=smpl_select(find(Subjpheno==1),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxTm_S=smpl_select(find(Subjpheno==0),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_A=smpl_select(find(Subjpheno==1),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T

samplmtrx_A=[samplmtrxT_A,samplmtrxTm_A];
samplmtrx_S=[samplmtrxT_S,samplmtrxTm_S];

samplmtrx_A(min(samplmtrx_A')'==0,:)=[]; %eliminate missing chips
samplmtrx_S(min(samplmtrx_S')'==0,:)=[];
Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);


[h,pv]=ttest2((data(ginds,samplmtrx_A(:,1))-data(ginds,samplmtrx_A(:,2)))',(data(ginds,samplmtrx_S(:,1))-data(ginds,samplmtrx_S(:,2)))',q,'both','unequal');
pvmax=max(pvmax,pv);

end

pv=(pvmax).^(length(timesel)-1);  

    


    case 'DeltaTphenomax'
      

pvmin=ones(1,Ngenes);
for T=1:length(timesel)-1 %minumum must be greater than 1                 
samplmtrxTm_A=smpl_select(find(Subjpheno==1),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxTm_S=smpl_select(find(Subjpheno==0),timesel(T),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_A=smpl_select(find(Subjpheno==1),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T
samplmtrxT_S=smpl_select(find(Subjpheno==0),timesel(T+1),excelmat,imtrx,alltimes,diagnostic);%Identify samples from 0.8T and T

samplmtrx_A=[samplmtrxT_A,samplmtrxTm_A];
samplmtrx_S=[samplmtrxT_S,samplmtrxTm_S];

samplmtrx_A(min(samplmtrx_A')'==0,:)=[]; %eliminate missing chips
samplmtrx_S(min(samplmtrx_S')'==0,:)=[];
Nsamp_A=length(samplmtrx_A);
Nsamp_S=length(samplmtrx_S);


[h,pv]=ttest2((data(ginds,samplmtrx_A(:,1))-data(ginds,samplmtrx_A(:,2)))',(data(ginds,samplmtrx_S(:,1))-data(ginds,samplmtrx_S(:,2)))',q,'both','unequal');
pvmin=min(pvmin,pv);

end

pv=1-(1-pvmin).^(length(timesel)-1);  


end %switch
end
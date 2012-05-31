function [Err,h]=OriginalVSReconstructedAligned(dopt,F,X,M,Masks,Padding,Genes,Subjects,figOnOff,SubjIDs,Title,GeneNames,OnsetTimes,Lambdas,TimeVec,GenePvs)
S=length(Subjects);
n=size(Masks{1},1);

k=size(F,2);
if figOnOff
h=figure('Name',strcat('Original VS Rec',Title));
end
Nplot=1;
Err=0;
Xr=cell(S,1);
for i=1:S
    ni=size(find(Masks{i}==1),1);
    Fadj=zeros(n,k);
    FadjAux=zeros(n,k);
    for j=1:k
         [T]=CreateTranslationMatrix(dopt{j}(Subjects(i)),n);
         Fadj(:,j)= T*F(:,j);
    end   
    n_i=size(M{Subjects(i)},1);
    TrsMask=(Masks{i}==1);
    Fadj=Fadj(TrsMask,:)*diag(sqrt(Lambdas));
    Xr{i}=real(F(Masks{i}==1,:)*diag(sqrt(Lambdas))*((Fadj'*Fadj)\Fadj')*X{Subjects(i)});
    Xr2{i}=real(Fadj*((Fadj'*Fadj)\Fadj')*X{Subjects(i)});
    XAxis=TimeVec;
    XAxis=XAxis(TrsMask);
    XAxis=XAxis((1+Padding):(n_i-Padding));
    if figOnOff
        subplot(ceil(S/2),4,Nplot)
        plot(XAxis,M{Subjects(i)}((1+Padding):(n_i-Padding),Genes),'LineWidth',2); hold on;
        %plot(TimeVec(OnsetTimes(Subjects(i))),max(max(Xr{i}(1+Padding:(n_i-Padding),Genes))),'ro','LineWidth',3); hold on        
        xMin=min(XAxis);
        xMax=max(XAxis);
        yMin=min(min(M{Subjects(i)}((1+Padding):(n_i-Padding),Genes)));
        yMax=max(max(M{Subjects(i)}((1+Padding):(n_i-Padding),Genes)));
        axis([xMin,xMax,1.1*yMin,1.1*yMax]);
        title(strcat('O (S. ID:',num2str(SubjIDs(i)),')'))
        subplot(ceil(S/2),4,Nplot+1)
        plot(XAxis,Xr2{i}(1+Padding:(n_i-Padding),Genes),'LineWidth',2);hold on;
%         if i==1
%            legend([GeneNames,'- Reconstructed Aligned'],'Reconstructed') 
%         end
        %plot(TimeVec(OnsetTimes(Subjects(i))),max(max((Xr{i}(1+Padding:(n_i-Padding),Genes)))),'ro','LineWidth',3); hold on
        axis([xMin,xMax,1.1*yMin,1.1*yMax]);
    end
    Err=Err+1/(S)*norm(Xr{i}-M{Subjects(i)},'fro')^2/norm(M{Subjects(i)},'fro')^2;
    if figOnOff
        title(strcat('R'))%, Err:', num2str(norm(M{Subjects(i)}-Xr{i},'fro')^2/norm(M{Subjects(i)},'fro')^2)))
    end
    if nargin==16
        if i==1
            lgnstr=cell(length(GeneNames),1);
            for l=1:length(GeneNames)
                lgnstr{l}=strcat(GeneNames{l},'(',num2str(GenePvs(l),'%10.1e\n'),')');
            end
            legend(lgnstr);
        end
        
    else
        if i==1
             legend(GeneNames);
         end
    end
    Nplot=Nplot+2;
end
suptitle(strcat('Average Relative Residual Error',num2str(Err)))

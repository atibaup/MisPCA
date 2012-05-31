function [Err,h]=OriginalVSReconstructed(dopt,F,X,M,Masks,Padding,Genes,Subjects,figOnOff,SubjIDs,Title,GeneNames,OnsetTimes,Lambdas,TimeVec,GenePvs)
S=length(Subjects);
n=size(Masks{1},1);

k=size(F,2);

if length(dopt)==S
    doptOld=dopt;
    dopt=cell(k,1);
   for i=1:k
       dopt{i}=[];
       for s=1:S
            dopt{i}=[dopt{i},doptOld{s}(i)];
       end
   end
end

Nplot=1;
Err=0;
Xr=cell(S,1);

NplotsPerFigure=10;

for i=1:S

    Fadj=zeros(n,k);
    for j=1:k
         Fadj(:,j)= circshift(F(:,j),-dopt{j}(Subjects(i)));
    end   

    MaskIndx=find(Masks{i}(:,1)==1);
    PlotIndx=intersect(MaskIndx,(1+Padding):(n-Padding));
   
    Fadj=Fadj*diag(sqrt(Lambdas));
    A=((Fadj'*Fadj)\(Fadj(MaskIndx,:)'*X{Subjects(i)}(MaskIndx,:)));

    Xr{i}=Fadj*A;

    signF=sign(sum(A,2));
    XAxis=TimeVec;
    XAxis=XAxis(PlotIndx);

    if figOnOff
        if mod(i-1,NplotsPerFigure)==0
             h=figure('Name',strcat(' ',Title));
             Nplot=1;
        end
        subplot(NplotsPerFigure/2,4,Nplot)
       
        OriginalData=X{Subjects(i)}(PlotIndx,Genes);
        RecoveredData= Xr{i}(:,Genes);
        
        plot(XAxis,OriginalData); hold on;
        plot(TimeVec,max(max(abs(OriginalData)))/max(max(abs(Fadj)))*Fadj*diag(signF),'LineWidth',3);
        plot(OnsetTimes(Subjects(i)),max(max(abs(OriginalData))),'ro','LineWidth',3); hold on
        
        title(strcat('Original (Subject ID:',num2str(SubjIDs(i)),')'))
        xMin=min(XAxis);
        xMax=max(XAxis);
        yMin=min(min(OriginalData));
        yMax=max(max(OriginalData));
        axis([xMin,xMax,1.1*yMin,1.1*yMax]);
        
        
        subplot(NplotsPerFigure/2,4,Nplot+1)
        plot(TimeVec,RecoveredData);hold on;
        stem(OnsetTimes(Subjects(i)),max(max((RecoveredData))),'ro'); hold on
        stem(TimeVec(mod(dopt{j}(Subjects(i)),n)+1),max(max((RecoveredData))),'bo'); hold on
        axis([xMin,xMax,1.1*yMin,1.1*yMax]);
    end
    Err=Err+1/(S)*norm(Xr{i}-M{Subjects(i)},'fro')^2/norm(M{Subjects(i)},'fro')^2;
    if figOnOff
        title(strcat('Reconstructed, Err:', num2str(norm(M{Subjects(i)}-Xr{i},'fro')^2/norm(M{Subjects(i)},'fro')^2)))
        
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
    end
    Nplot=Nplot+2;
end
if figOnOff
   strcat('Average Relative Residual Error =',num2str(Err))
end

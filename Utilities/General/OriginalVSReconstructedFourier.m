function [Err,h]=OriginalVSReconstructedFourier(dopt,F,X,For,dor,Ft,figOnOff,Title,GeneNames)
S=length(X);
n=size(X{1},1);
p=size(X{1},2);
k=size(F,2);
if figOnOff
h=figure('Name',strcat('Original VS Rec',Title));
end
DFT=1/sqrt(n)*dftmtx(n);
DFTinv=conj(DFT);
Nplot=1;
Err=0;
for i=1:S
    %     Fadj=zeros(n,k);
    %     for j=1:k
    %         D=diag(exp(-1i*2*pi/n*[0:n-1]*dopt(i)));
    %         Fadj(:,j)=D*DFT*F(:,j);
    %     end
    % Fadj=(DFTinv*Fadj);
    %Xr{i}=(Fadj*((Fadj'*Fadj)\(Fadj'*X{Subjects(i)})));
    if figOnOff
        Dor=diag(exp(-1i*2*pi/n*[0:n-1]*dor{i}(1)));
        Hor=DFTinv*Dor*DFT*For(:,1);
        Dopt=diag(exp(-1i*2*pi/n*[0:n-1]*dopt{1}(i)));
        Hopt=DFTinv*Dopt*F(:,1);
        
        subplot(S,4,Nplot)
        plot(abs(Hor)); hold on
        %plot(OnsetTimes(Subjects(i)),max(max((Xr{i}(1+Padding:(n-Padding),Genes)))),'ro','LineWidth',3); hold on
        title(strcat('Original (mod)'))
        axis tight;
        
        subplot(S,4,Nplot+1)
        plot(abs(Hopt));hold on;
        %plot(OnsetTimes(Subjects(i)),max(max((Xr{i}(1+Padding:(n-Padding),Genes)))),'ro','LineWidth',3); hold on
        title(strcat('Recovered (mod)'))
        axis tight;
        
        subplot(S,4,Nplot+2)
        plot(angle(Hor)); hold on
        %plot(OnsetTimes(Subjects(i)),max(max((Xr{i}(1+Padding:(n-Padding),Genes)))),'ro','LineWidth',3); hold on
        title(strcat('Original (phase)'))
        axis tight;
        
        subplot(S,4,Nplot+3)
        plot(angle(Hopt));hold on;
        %plot(OnsetTimes(Subjects(i)),max(max((Xr{i}(1+Padding:(n-Padding),Genes)))),'ro','LineWidth',3); hold on
        title(strcat('Recovered (mod)'))
        axis tight;
        
    end
    %     if figOnOff
    %         title(strcat('Reconstructed (mod), Err:', num2str(norm(M{Subjects(i)}-Xr{i},'fro'))))
    %     end
    %     if figOnOff
    %         title(strcat('Reconstructed (phase), Err:', num2str(norm(M{Subjects(i)}-Xr{i},'fro'))))
    %     end
    if nargin==11
        if i==1
            legend(GeneNames);
        end
    end
    Err=Err+1/(S)*norm(Hopt-Hor,'fro')^2;
    
    Nplot=Nplot+4;
end

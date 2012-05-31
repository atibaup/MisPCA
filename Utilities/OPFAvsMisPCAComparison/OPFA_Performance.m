function [Fr,F,FactorsDistance,ErrM,ErrX]=OPFA_Performance(X,D,Output,F,d,FiguresOnOff,titlestr)
n_F=size(F,1);
n=size(X{1},1);
S=length(X);
f_o=size(Output.F,2);
epsilon=n_F-n;
epsilon1=round(epsilon/2);

ErrM=0;
ErrX=0;

if FiguresOnOff
    h=figure('Name',strcat('Reconstructed vs Original (trajectories)',titlestr));
    h2=figure('Name',strcat('Reconstructed vs Original (factors)',titlestr));
end
kS=1;
Fr=cell(S,1);
Fo=cell(S,1);
FactorsDistance=zeros(S,1);

for i=1:S
    Mr=[];
    for j=1:f_o
        AlignedFact=circshift(Output.F(:,j),-Output.dhat{i}(j));
        Mr=[Mr AlignedFact(epsilon1+1:epsilon1+n)];
    end
    M=[];
    for j=1:f_o
        AlignedFact=circshift(F(:,j),-d{i}(j));
        M=[M AlignedFact(epsilon1+1:epsilon1+n)];
    end
    
    Fr{i}=Mr;
    Fo{i}=M;

    FactorsDistance(i)=FactorsSubspaceDistance(Fr{i},Fo{i});
    
    Rec=Mr*Output.A_o{i};
    
    if FiguresOnOff
        figure(h2)
        subplot(S,3,kS)
        plot(D{i})
        if kS==1; title('Original'); end;
        subplot(S,3,kS+1)
        plot(Rec)
        if kS==1; title('Reconstructed'); end;
        subplot(S,3,kS+2)
        plot(D{i}-Rec)
        if kS==1; title('Residual'); end;
        kS=kS+3;
        figure(h)
        subplot(2,1,1)
        plot(Fr{i});hold on;
        title('Reconstructed')
        subplot(2,1,2)
        plot(Fo{i});hold on;
        title('Original')
    end
    ErrM=ErrM+1/S*norm(D{i}-Rec,'fro')/norm(D{i},'fro');
    ErrX=ErrX+1/S*norm(X{i}-Rec,'fro')/norm(X{i},'fro');
end
if FiguresOnOff
    if exist('h2','var')
        set(h2,'Name',strcat('Reconstructed vs Original (trajectories)',titlestr,'-',['Average MSE =',num2str(ErrM)]));
        figure(h2);suptitle(strcat('Reconstructed vs Original (trajectories)',titlestr,'-',['Average MSE =',num2str(ErrM)]));
    end
    if exist('h','var')
        set(h,'Name',strcat('Reconstructed vs Original (factors)',titlestr,'-',['Average DTF =',num2str(mean(FactorsDistance))]));
        figure(h);suptitle(strcat('Reconstructed vs Original (factors)',titlestr,'-',['Average DTF =',num2str(mean(FactorsDistance))]));
    end
end
FactorsDistance=mean(FactorsDistance);

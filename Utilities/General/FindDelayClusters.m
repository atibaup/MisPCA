function [Clusters]=FindDelayClusters(S,n_f,p,X,k,FigOnOff)

VecData=zeros(S,n_f*p);
for i=1:S
    VecData(i,:)=abs(X{i}(:));
end

T=clusterdata(VecData,'maxclust',k,'distance','correlation','linkage','complete');

[Vals,indxsorted]=sort(T,'ascend');
if FigOnOff
    figure
    subplot(1,2,1)
    imagesc(VecData(indxsorted,:))
    subplot(1,2,2)
    imagesc(Vals)
end
for i=1:k
    Indx=find(T==i);
    if FigOnOff
        figure
        for j=1:length(Indx)
            subplot(length(Indx),1,j)
            plot(X{Indx(j)});hold on;
%            plot(M{Indx(j)},'LineWidth',3);hold on;
        end
    end
    Clusters{i}=Indx;
end

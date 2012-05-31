function [Xaligned,Aaligned]=AlignedRepresentation(X,Masks,Results)

S=length(X);


d=Results.d;
F=Results.F;
Lambdas=Results.Lambdas;
f=size(F,2);
n=size(F,1);
dRef=round(mean(d{1}));

for i=1:S
    for j=1:f
       dInd{i}(j)= d{i}(j);
    end
end
for i=1:S
    Fadj=zeros(n,f);
    Fal=zeros(n,f);
    for j=1:f

        Fadj(:,j)=circshift(F(:,j),dInd{i}(j));
        Fal(:,j)=circshift(F(:,j),dInd{i}(j)-dInd{i}(1)+dRef);
    end
    Fadj=Fadj(Masks{i}(:,1)==1,:)*diag(sqrt(Lambdas));
    Aaligned{i}=Fadj\X{i}(Masks{i}(:,1)==1,:);
    Xaligned{i}=Fal*Aaligned{i};
end
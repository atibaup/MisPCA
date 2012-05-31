function [newX,NewMasks,newy,commonTimeVec]=combine_PHD_datasets(X,Masks,TimeVecOr,y,Padding)
Nds=length(X);
commonTimeVec=[];
for i=1:Nds
   commonTimeVec=union(commonTimeVec,round(TimeVecOr{i})); 
end
commonTimeVec=commonTimeVec
newN=length(commonTimeVec);
k=1;
newy=[];
for i=1:Nds
   [C,IA,IB]=intersect(commonTimeVec,round(TimeVecOr{i}));
   
    for j=1:length(Masks{i})
        NewMasks{k}=zeros(newN,size(Masks{i}{j},2));
        newX{k}=zeros(newN,size(Masks{i}{j},2));
        NewMasks{k}(IA,:)=Masks{i}{j};

        newX{k}(IA,:)=X{i}{j};
        
        if Padding>0
            p=size(NewMasks{k},2);
            newX{k}=[repmat(newX{k}(1,:),Padding+newN-size(newX{k},1),1); newX{k}; repmat(newX{k}(end,:),Padding+newN-size(newX{k},1),1)];
            NewMasks{k}=[ones(Padding,p);NewMasks{k};ones(Padding,p)];
        end
        k=k+1;
    end
    newy=[newy;y{i}(:)];
end
    
if Padding>0
   commonTimeVec=[(-Padding:-1)+min(commonTimeVec),commonTimeVec,max(commonTimeVec)+(1:Padding)];
end
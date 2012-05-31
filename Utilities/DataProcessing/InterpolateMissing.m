function [X,Mask]=InterpolateMissing(X,Mask)

S=length(X);


for i=1:S
    fprintf('\n Interpolating subject %d/%d',i,S);
    p=size(X{i},2);
    n=size(X{i},1);
    x=find(Mask{i}(:,1)==1);
    if length(x)<n
        xi=find(Mask{i}(:,1)==0);
        Xinterp{i}=zeros(n,p);
        for g=1:p
            Y=X{i}(x,g);
            Xinterp{i}(x,g)=Y;
            Xinterp{i}(xi,g)=interp1(x,Y,xi);
        end
        Mask{i}=ones(n,p);
        X{i}=Xinterp{i};
    end
end
fprintf('\n');
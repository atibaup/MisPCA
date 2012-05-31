function [c]=constructC(Sigma,H,Set,GridSize,dmax)
n=size(Sigma{1},1);

Grid=0:GridSize:dmax;
nGrid=length(Grid);
S=length(Sigma);
c=zeros(nGrid*S,1);
for i=1:S
    for j=1:nGrid
        if j-1>= Set(i,1) && j-1 <= Set(i,2)
            T=CreateTranslationMatrix(Grid(j),n);
            c((i-1)*nGrid+j)=trace(H*T'*Sigma{i}*T);
        end
    end
end

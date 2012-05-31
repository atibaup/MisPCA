function [S]=GenerateOPPCAData(F,d,NS,p,n)
f=length(d{1});
for i=1:NS
   X=randn(p,n);
   Fal=zeros(p,f);
    for j=1:f
            [T]=CreateTranslationMatrix(d{i}(j),p);
            Fal(:,j)=T*F(:,j);
    end
   S{i}=Fal*Fal'+1/n*X*X';
end
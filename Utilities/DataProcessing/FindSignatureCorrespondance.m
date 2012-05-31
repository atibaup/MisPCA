function [Corr]=FindSignatureCorrespondance(Signature1,Signature2)

n=length(Signature1);
Left=1:n;
for i=1:n
   Cost=zeros(n,1); 
   for j=1:length(Left)
   
       Cost(j)=max(xcorr(Signature1{i}.AverageSign,Signature2{Left(j)}.AverageSign));
   end
   [maxC,Indx]=max(Cost);
   Corr(i)=Left(Indx(1));
   Left=setdiff(Left,Left(Indx(1)));
end
   
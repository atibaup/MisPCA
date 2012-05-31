function [Grid]=CreateNdimGrid(Set)
f=size(Set,1);
d(1)=1;
r(f)=1;
for i=2:f 
    d(i)=d(i-1)*(max(Set(i-1,:))-min(Set(i-1,:))+1);
    r(f-i+1)=r(f-i+2)*(max(Set(i,:))-min(Set(i,:))+1);
end

Ngrid=d(f)*(max(Set(f,:))-min(Set(f,:))+1);

Grid=zeros(Ngrid,f);
for i=0:f-1
    for j=1:d(f-i)
%         length((j-1)*d(f-i)+1:(j*d(f-i)))
%         length(min(Set(:,f-i)):max(Set(:,f-i)))
        Aux=RepVector(min(Set(f-i,:)):max(Set(f-i,:)),r(f-i));%min(Set(:,f-i)):max(Set(:,f-i));
        Grid((j-1)*length(Aux)+1:(j*length(Aux)),f-i)=Aux;
    end
    
end

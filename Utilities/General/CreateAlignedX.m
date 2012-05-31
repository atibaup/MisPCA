function [X_a]=CreateAlignedX(X,d)

S=length(d);
n=size(X{1},2);
X_a=zeros(n);
for i=1:S
    [T]=CreateTranslationMatrix(d(i),n);
    X_a=X_a+T'*X{i}*T;
end
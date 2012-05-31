function [LB,UB]=ToeplitzEigMaxBounds(h)

h=h(:)';
n=length(h);
hflip=fliplr(h(1:end));
has=(h+hflip);
symmh=[h,[0 hflip(1:end-1)]];

W=dftmtx(2*n);
predict=(W*symmh');
predict=sort(predict,'descend');
predictEven=predict(1:2:end);
predictOdd=predict(2:2:end);
UB=abs(1/2*(predictEven(1)+predictOdd(1)));
LB=abs(1/2*(predictEven(1)+predictOdd(end)));
function [dMSE,ProjMSE,MSE,bias]=ProjAndDelayMSE(F_o,d_o,Fhat,dhat)

n=size(F_o,1);
S=length(d_o);
if iscell(d_o)
    for i=1:S
        dvec(i)=d_o{i};
    end
else
    dvec=d_o;
end

%dMSE=1/S*norm((dvec-min(dvec))-(dhat-min(dhat)))^2;
dMSE=0;
ProjMSE=0;
MSE=0;
bias=0;
%figure
for i=1:S
    h=CreateTranslationMatrix(dvec(i),n)*F_o;
    hhat=CreateTranslationMatrix(dhat(i),n)*Fhat;
    ProjMSE=ProjMSE+1/S*abs((1/norm(h,2)*h'*1/norm(hhat,2)*hhat));
    % plot(h,'-ob');hold on;
    % plot(hhat,'-xr');hold on;
    MSE=MSE+1/S*min(norm(h-hhat)^2,norm(h+hhat)^2);
    bias=bias + +1/S*min(mean(h-hhat),mean(h+hhat));
end


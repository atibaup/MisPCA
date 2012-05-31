function [y]=smoothSig(x,order,type)

if strcmp(type,'triang')
Coeffs=triang(order);
elseif strcmp(type,'rect')
    Coeffs=rectwin(order);
elseif strcmp(type,'exp')
    Coeffs=exp(-1/10*1:order);
    Coeffs=2/norm(Coeffs)*Coeffs;
    Coeffs=Coeffs(:);
end


Coeffs=[Coeffs;zeros(length(x)-length(Coeffs),1)];

y=cconv(x(:),Coeffs(:));
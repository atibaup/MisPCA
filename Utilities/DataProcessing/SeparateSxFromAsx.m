function [A_oAsx,A_noAsx,A_oSx,A_noSx,XSx,XAsx,kSx,kAsx]=SeparateSxFromAsx(A_o,A_no,X,y,Order,S)
kAsx=1;
kSx=1;
clear A_oAsx A_noAsx XSx A_oSx A_noSx XAsx
A_oAsx=cell(length(find(y==1)),1);
A_noAsx=cell(length(find(y==1)),1);
XAsx=cell(length(find(y==1)),1);
A_oSx=cell(length(find(y==0)),1);
A_noSx=cell(length(find(y==0)),1);
XSx=cell(length(find(y==0)),1);
for i=1:S
    if y(Order(i))==1
        A_oAsx{kAsx} =A_o{Order(i)}; 
        A_noAsx{kAsx} =A_no{Order(i)};
        XAsx{kAsx}=X{Order(i)};
        kAsx=kAsx+1;
    else
        A_oSx{kSx} =A_o{Order(i)}; 
        A_noSx{kSx} =A_no{Order(i)};
        XSx{kSx}=X{Order(i)};
        kSx=kSx+1;
    end
end
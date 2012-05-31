function [c,a]=constructC_v2(Sigma,Masks,H,Set,GridSize)
n=size(Masks{1},1);
S=length(Sigma);
c=zeros(n*S,1);
%tic
Niter=0;
a=zeros(n*S,1);
for i=1:S
    %ShowGrid= Set(i,1)+1:GridSize:Set(i,2)+1
        maxSoFar=0;    
    for j= Set(i,1)+1:GridSize:Set(i,2)+1

        %       if j-1>= Set(i,1) && j-1 <= Set(i,2)
        %T=CreateTranslationMatrix(j-1,n);
        %             Aux=zeros(n);
        %             Aux(Masks{i}==1,Masks{i}==1)=Sigma{i};
        %             Aux=T'*Aux*T;
        %             TrsMask=T*Masks{i};
        %             TrsMask=(TrsMask==1);
        %             Aux=Aux(TrsMask,TrsMask);
        %             c((i-1)*n+j)=trace(H(TrsMask,TrsMask)*Aux)
        TransH=circshift((circshift(H',-(j-1)))',-(j-1));
        %TransH=T*H*T';
        %Aux=Sigma{i};
        c((i-1)*n+j)=trace(TransH(Masks{i}==1,Masks{i}==1)*Sigma{i});
        %       end
        if c((i-1)*n+j)>=maxSoFar
%             'enter'
%             Indx=(i-1)*n+j
%             currmax=c((i-1)*n+j)
%             pastmax=maxSoFar
           a((i-1)*n+1:i*n)=zeros(n,1);
           a((i-1)*n+j)=1;
           maxSoFar=c((i-1)*n+j);
        end
        Niter=Niter+1;
    end
end
Niter=Niter
%Timeloop=toc
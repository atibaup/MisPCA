function [Err,Merr]=PredictionErr(X,Masks,OmegaT,Xhat,FiguresOnOff)

S=length(X);

Err=zeros(S,1);
if FiguresOnOff;
figure;
end

k=1;

for i=1:S

    MaskIndx=(Masks{i}(:,1)==1);
    Xor=X{i}(MaskIndx,:);
    Xrec=Xhat{i}(MaskIndx,:);   
    
    TestIndx=(OmegaT{i}(MaskIndx==1,:)==1);
    
    Err(i)=1/length(TestIndx)*norm(Xor(TestIndx)-Xrec(TestIndx),'fro')^2;
    
    if FiguresOnOff;       
        subplot(S,2,k)
        plot(Xor)
        subplot(S,2,k+1)
        plot(Xrec)
        title(['Test error=',num2str(Err(i))])
        k=k+2;
    end

end
Merr=mean(Err);
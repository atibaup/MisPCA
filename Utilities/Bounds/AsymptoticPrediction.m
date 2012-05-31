function [PredictedLambda,PredictedCorr,B]=AsymptoticPrediction(H,d,s_d,sigma,SNR,FigOnOff)
n=length(d);
F=size(H,2);
p=size(H,1);
Npoints=length(SNR);

SNRdB=SNR;
SNR=10.^(SNR/10);

[B,PT,PT2,dummy1,dummy2,lambda]=PhaseTransitionBounds(H,d,s_d,sigma);
 

PredictedCorr=zeros(Npoints,F);
PredictedLambda=zeros(Npoints,F);
c=p/n;

for f=1:F
    
    gammaf=sqrt(c)/lambda(f);
    PredictedLambda(SNR<=gammaf,f)=(1+sqrt(c))^2;
    PredictedCorr(SNR<=gammaf,f)=0;
    
    SNRpl=SNR(SNR>gammaf);
    PredictedLambda(SNR>gammaf,f)=(SNRpl*lambda(f)+1).*(1+c./(SNRpl*lambda(f)));
    PredictedCorr(SNR>gammaf,f)=((SNRpl*lambda(f)).^2-c)./((SNRpl*lambda(f)).^2+c*SNRpl*lambda(f));
    
end


if FigOnOff
    figure;
    subplot(2,1,1);
    plot(SNRdB,PredictedLambda,'-o'); hold on;
    %stem(PT,max(PredictedLambda),'-x'); hold on;
    xlabel('SNR (dB)')
    ylabel('\lambda_f(S)')
    
    subplot(2,1,2);
    plot(SNRdB,PredictedCorr,'-o'); hold on;
    %stem(PT,max(PredictedCorr),'-x'); hold on;
    xlabel('SNR (db)')
    ylabel('|<h_f ,v_f(S) >|^2')
end

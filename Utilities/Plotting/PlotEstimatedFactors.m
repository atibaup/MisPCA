function []=PlotEstimatedFactors(F_MisPCA,dopt_MisPCA,Lambdas,Options,OnsetTimes,ChosenFactors,TimeVec)
h2=figure('Name','Estimated Factors')
colors=hsv(length(ChosenFactors))
markers=['+o.*xs^v']
for i=1:length(ChosenFactors)
    subplot(ceil(length(ChosenFactors)/2)+1,2,i)
%    plot(F_MisPCA(Options.Padding+1:(Options.n-Options.Padding),i),'Color',colors(i,:),'LineWidth',2); hold on;
    for j=1:Options.S
        T=CreateTranslationMatrix(dopt_MisPCA{ChosenFactors(i)}(j),Options.n);
        Fadj=T*F_MisPCA(:,ChosenFactors(i));
        if j==1
        plot(Fadj(Options.Padding+1:(Options.n-Options.Padding)),'Color','r','LineWidth',3); hold on;
        else
            plot(Fadj(Options.Padding+1:(Options.n-Options.Padding)),'Color','r'); hold on;
        end
        stem(OnsetTimes(j),max(Fadj(Options.Padding+1:(Options.n-Options.Padding))),'x','Color','k');hold on;
    end
    title(strcat('Factor ',num2str(ChosenFactors(i))))
end
    suptitle('Estimated Factors after translation')
if mod(length(ChosenFactors),2)==0
subplot(ceil(length(ChosenFactors)/2)+1,2,i+[1,2])
else
    subplot(ceil(length(ChosenFactors)/2)+1,2,i+1+[1,2])
end
plot(Lambdas(ChosenFactors),'-o','LineWidth',2)
title('Eigenvals')
% for i=1:length(ChosenFactors)
%     h3=figure;
%     for j=1:Options.S
%         subplot(ceil(Options.S/2),2,j)
%         T=CreateTranslationMatrix(dopt_MisPCA{ChosenFactors(i)}(j),Options.n);
%         Fadj=T*F_MisPCA(:,ChosenFactors(i));
%         plot(TimeVec(Options.Padding+1:(Options.n-Options.Padding)),Fadj(Options.Padding+1:(Options.n-Options.Padding)),'Color','r'); hold on;
%         stem(TimeVec(OnsetTimes(j)),max(Fadj(Options.Padding+1:(Options.n-Options.Padding))),'x','Color','k');hold on;
%     end
%     suptitle(strcat('Factor profiles for factor:',num2str(ChosenFactors(i))));
% end
function [AdjustedDelays]=PlotEstimatedDelays(X,OnsetTimes,dopt_MisPCA,RefFactor,Options,ChosenFactors)
PeakTimes=OnsetTimes;
h3=figure('Name','Estimated Delays');
markers=['+o.*xs^v'];
Nfact=length(ChosenFactors);
colors=cool(Nfact);
[dummy, SorteddoptIndx]=sort(dopt_MisPCA{RefFactor});
AdjustedDelays=cell(Nfact,1);
Lgndstr=cell(Nfact,1);
for i=1:Nfact
    AdjustedDelays{i}=mod(Options.Padding+dopt_MisPCA{ChosenFactors(i)}(:)-max(dopt_MisPCA{ChosenFactors(i)}(:)),Options.n);
   % plot(1:length(X), mod(PeakTimes(:)-Options.Padding+dopt_MisPCA{i}(SorteddoptIndx),Options.n),'LineStyle','-',...
    %    'Marker',markers(i),'Color',colors(i,:)); hold on;
    plot(1:length(X), mod(dopt_MisPCA{ChosenFactors(i)}(SorteddoptIndx),Options.n),'LineStyle','-',...
        'Marker',markers(i),'Color',colors(i,:),'LineWidth',2); hold on;
    Lgndstr{i}=['Peak Factor \#',num2str(i)];
end
plot(1:length(X), OnsetTimes(SorteddoptIndx),'LineStyle','-' ,'Marker',markers(i),'Color',colors(i,:),'LineWidth',3); hold on;
Lgndstr{i+1}='Subject Onset Times';
legend(Lgndstr)

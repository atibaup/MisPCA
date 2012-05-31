function [delays,dAbsolute]=MatchDelays(F,M_no,dhats,SxorAsx,OnsetTimes,ChangePoint,RefFactor,TimeVec,FigOnOff)
Mi=[];
S=length(dhats);
f_o=size(F,2);
n_F=size(F,1);
f_no=size(M_no,2);
n=size(M_no,1);
epsilon=n_F-n;
epsilon1=round(epsilon/2);
epsilon2=epsilon-epsilon1;
for i=1:S
    for j=1:f_o
        Ui=CreateTranslationMatrix(dhats{i}(j),n_F);
        AlignedFact=Ui*F(:,j);
        Mi=[Mi ,AlignedFact(epsilon1+1:epsilon1+n)];
    end
end
if SxorAsx==1
    %[m FirstPeak]=max(abs(Mi(:,1)))
    FirstPeak=ChangePoint;
else
    %[m FirstPeak]=max(abs(Mi(:,2)))
    FirstPeak=ChangePoint;
end
diff=FirstPeak-dhats{1}(RefFactor);
del=dhats{1}(1);
%m=m
delays=[];
ActualDelays=[];
t=1;
ticks(1)={''};
for i=1:S
    for j=1:f_o
        delays=[delays TimeVec(mod(dhats{i}(j)+diff,n_F)+1)];
        dAbsolute{i}(j)=mod(dhats{i}(j)+diff,n);
         ticks(t)={strcat('S',num2str(i),'-F',num2str(j))};
        t=t+1;
        if iscell(OnsetTimes)==0
            ActualDelays=[ActualDelays TimeVec(OnsetTimes(i))];
        else
            ActualDelays=[ActualDelays OnsetTimes{i}(j)];
        end
    end
end


% figure('Name','Matched Delays and Factors')
% subplot(2,f_o+f_no,[1:f_o+f_no])
% imagesc(1:f_o*S,TimeVec,Mi); hold on;
% h=plot(delays,'-*r','LineWidth',4);hold on;
% plot(ActualDelays,'-*g','LineWidth',4)
% set(gca,'XTick',1:(t-1), 'XTickLabel',ticks,'FontSize',12,'FontWeight','bold')
% ylabel('Time (hours)','FontSize',12,'FontWeight','bold')
% xlabel('Subject - Factor','FontSize',12,'FontWeight','bold')
%    axis([0 S*f_o min(TimeVec) max(TimeVec)])
% for i=1:f_o
%    subplot(2,f_o+f_no,(f_o+f_no)+i);
%    plot(F(:,i),'-*r','LineWidth',4)
%    xlabel('Time','FontSize',12,'FontWeight','bold')
%    ylabel('Expression Value','FontSize',12,'FontWeight','bold')
%    title(strcat('Factor:',num2str(i)),'FontSize',12,'FontWeight','bold')
%    axis tight
% end
% for i=1:f_no
%    subplot(2,f_o+f_no,2*f_o+f_no+i);
%    plot(M_no(:,i),'-*r','LineWidth',4)
%    xlabel('Time','FontSize',12,'FontWeight','bold')
%    ylabel('Expression Value','FontSize',12,'FontWeight','bold')
%    title(strcat('NO Factor:',num2str(i)),'FontSize',12,'FontWeight','bold')
%       axis tight
% end
if FigOnOff==1
    figure('Name','Matched Delays')
    subplot(2,f_o+f_no,[1:f_o+f_no])
    h=plot(delays,'-sr','LineWidth',4,'MarkerSize',10);hold on;
    plot(ActualDelays,'-og','LineWidth',4,'MarkerSize',10)
    set(gca,'XTick',1:(t-1), 'XTickLabel',ticks,'FontSize',12,'FontWeight','bold')
    ylabel('Time (hours)','FontSize',12,'FontWeight','bold')
    xlabel('Subject - Factor','FontSize',12,'FontWeight','bold')
    legend('Estimated Factor Delays','Peak Symptom Time')
    axis([0 S*f_o min(TimeVec) max(TimeVec)])
    for i=1:f_o
        subplot(2,f_o+f_no,f_o+f_no+i);
        plot(F(:,i),'-*r','LineWidth',4)
        xlabel('Time','FontSize',12,'FontWeight','bold')
        ylabel('Expression Value','FontSize',12,'FontWeight','bold')
        title(strcat('Factor:',num2str(i)),'FontSize',12,'FontWeight','bold')
        axis tight
    end
    
    for i=1:f_no
        subplot(2,f_o+f_no,2*f_o+f_no+i);
        plot(M_no(:,i),'-*r','LineWidth',4)
        xlabel('Time','FontSize',12,'FontWeight','bold')
        ylabel('Expression Value','FontSize',12,'FontWeight','bold')
        title(strcat('NO Factor:',num2str(i)),'FontSize',12,'FontWeight','bold')
        axis tight
    end
end
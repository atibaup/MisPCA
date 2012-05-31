function PlotMisPCAwCVResults(Results)

figure
subplot(4,Results.f,[1:Results.f])
imagesc(Results.CVErr); hold on;
plot(Results.iopt,Results.kopt,'or','MarkerSize',8);
xlabel('\beta')
ylabel('f')
title('CV error results')

colors=hsv(Results.f);
for k=1:Results.f
    subplot(4,Results.f,Results.f+k)
    
    plot(circshift(Results.F(:,k),Results.d{k}(1)),'color',colors(k,:),'Linewidth',2);
    %f=length(Results.d);
    % legendstr=cell(f,1);
    % for i=1:f
    %    legendstr{i}=['Factor ',num2str(i)];
    % end
    % legend(legendstr)
    xlabel('Time')
    ylabel('Magnitude')
    axis tight;
    title(['MisPCA Factor ',num2str(k)])
end
subplot(4,Results.f,2*Results.f+ [1:Results.f])
marker=['xos*v^xos*v^'];
for i=1:Results.f
   plot(Results.d{i},'color',colors(i,:),'Linewidth',2,'Marker',marker(i)); hold on; 
end
xlabel('Subject')
ylabel('Time')
title('MisPCA delays')

subplot(4,Results.f,3*Results.f+ [1:Results.f])
plot(Results.Lambdas)
xlabel('Factor')
title('MisPCA eigvalss')

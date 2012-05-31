function [A_o,U,F,X,M,dini,SNRest,Ft,Masks,dmax]=GenerateFMisPCADataForOPFA(p,n,n_f,S,f_o,SNR,VarDelays,FigOnOff)

epsilon=n_f-n;
epsilon1=round(epsilon/2);


HalfIndx=round(n_f*.3);
Rest=round(1/2*(n_f-HalfIndx));
HalfIndx=n_f-2*Rest;
MotifDictionary=zeros(n_f,f_o);

MotifDictionary(:,1)=[ones(Rest,1);zeros(HalfIndx,1);zeros(Rest,1);]; % upregulation
MotifDictionary(:,2)=circshift(flipud([1.5*ones(round(HalfIndx/4),1);1.5*zeros(n_f-round(HalfIndx/4),1)]),-round(HalfIndx)); % downregulation
%MotifDictionary(:,3)=flipud([3*zeros(2*Rest-1,1);3*ones(1,1)+3*triang(1);3*zeros(HalfIndx,1)]);% upregulation peak
for i=3:f_o
    if rand>.5
         MotifDictionary(:,i)=smoothSig(circshift(sort(abs(randn(n_f,1))),round(rand*n_f))',round(3+rand*10),'triang');
    else
        MotifDictionary(:,i)=smoothSig(fliplr(circshift(sort(abs(randn(n_f,1))),round(rand*n_f))'),round(3+rand*10),'triang');
    end
end

for i=1:f_o
    MotifDictionary(:,i)=smoothSig(MotifDictionary(:,i)',6,'triang');
end


%% Generate Factors
F=zeros(n_f,f_o);
F_=zeros(n_f*f_o,f_o);
for i=1:f_o
    F(:,i)=abs(MotifDictionary(:,i)+ 1/100*randn(n_f,1));
    F(:,i)=1/norm(F(:,i))*F(:,i);
    F_((i-1)*n_f+1:i*n_f,i)=F(:,i);
end

if FigOnOff
    figure('Name','Data model')
    cols=hsv(f_o);
    %figure('Name','Factors')
    legstr=cell(f_o,1);
    for i=1:f_o
        plot(F(:,i),'Linewidth',3,'Color',cols(i,:));hold on;
        xlabel('Time','FontWeight','bold','FontSize',12)
        ylabel('Magnitude','FontWeight','bold','FontSize',12)
        legstr{i}=['Factor',blanks(1),num2str(i)];
    end
    legend(legstr)
end

%% Generate subject's data
SNRest=0;
dmax=0;
A_o=cell(S,1);
dini=cell(S,1);
for w=1:S

    dini{w}=(floor(rand*(round(sqrt(12*VarDelays+1))))+1)*ones(f_o,1);
    dmax=max(dmax,max(dini{w}));
    U{w}=[];
    t=0:n-1;
    AlignedFactors=zeros(n_f,f_o);
    for q=1:f_o
        if dini{w}(q)>(n_f-1)
            dini{w}(q)=n_f-1;
        end
        AlignedFactors(:,q)=circshift(F(:,q),-dini{w}(q));
        U{w}=[U{w} sparse(circshift(eye(n_f),-dini{w}(q)))];
    end
    Ft{w}=AlignedFactors;
    A_o{w}=abs(randn(f_o,p));    
    
    Data=AlignedFactors(epsilon1+1:epsilon1+n,:)*A_o{w};
    Noise=1/sqrt(n*p)*sqrt(norm(Data,'fro')^2/(10^(SNR/10)))*randn(n,p);
    X{w}=Data+Noise;
    M{w}=Data;
    SNRest=SNRest+norm(Data,'fro')/norm(Noise,'fro');
end
SNRest=20*log10(1/S*SNRest);
if FigOnOff
    h4=figure('Name','Generated Data');
    for i=1:S
        figure(h4)
        subplot(2,ceil(S/2),i)
        plot(X{i})
        title(['Subject ',num2str(i),'/',num2str(S)]);
        xlabel('Time')
        ylabel('Magnitude')
        %figure
        %plot(X{i})
    end
end
Masks=cell(S,1);
for i=1:S
    Masks{i}=ones(n_f,p);
end

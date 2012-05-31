function [X,Sigma,SigmaAv,M,A,F,dini,SNRest,Ft,Masks,dmax]=GenerateMisPCAData(p,n,n_f,S,f_o,SNR,dmax,Spacing,FigOnOff,MissingTPRatio)

epsilon=n_f-n;
epsilon1=round(epsilon/2);


HalfIndx=round(n_f*.8);
Rest=round(1/2*(n_f-HalfIndx));
HalfIndx=n_f-2*Rest;
MotifDictionary=zeros(n_f,f_o);

MotifDictionary(:,1)=[ones(Rest,1);zeros(HalfIndx,1);zeros(Rest,1);]; % upregulation
MotifDictionary(:,2)=circshift(flipud([1.5*ones(round(HalfIndx/4),1);1.5*zeros(n_f-round(HalfIndx/4),1)]),-round(HalfIndx)); % downregulation
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

MotifDictionary=MotifDictionary-repmat(mean(MotifDictionary,1),n_f,1);
MotifDictionary=MotifDictionary*diag(sqrt(diag(MotifDictionary'*MotifDictionary)));
    
%% Generate motif-gene association
Pattern=zeros(p,f_o);
LeftIndx=1:p;
for i=1:f_o
    if i<f_o
         RandHalf=(1/2*length(LeftIndx));
     else
         RandHalf=length(LeftIndx);
    end
    RandHalf=round(RandHalf);
    CurrIndx=LeftIndx(1:RandHalf);
    LeftIndx=setdiff(LeftIndx,CurrIndx);
    Pattern(CurrIndx,i)=ones(length(CurrIndx),1);
end

%% Generate Factors
F=zeros(n_f,f_o);
F_=zeros(n_f*f_o,f_o);
for i=1:f_o
    F(:,i)=MotifDictionary(:,i)+ 1/100*randn(n_f,1);
    F(:,i)=1/norm(F(:,i))*F(:,i);
    F_((i-1)*n_f+1:i*n_f,i)=F(:,i);
end

if FigOnOff
    figure('Name','Data model')
    subplot(2,1,1)
    imagesc(Pattern)
    xlabel('Scores')
    ylabel('Factor')
    title('Scores sparsity pattern')
    subplot(2,1,2)
    cols=hsv(f_o);
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
Masks=cell(S,1);
dini=cell(S,1);
Ft=cell(S,1);
A=cell(S,1);
X=cell(S,1);
M=cell(S,1);
for w=1:S

    delays=unidrnd(round((dmax+1)/Spacing),1,f_o)-1;
    if dmax>0
        SpacedDelays=0:Spacing:dmax;
        dini{w}=SpacedDelays(delays+1);
    else
        dini{w}=0;
    end
    
    AlignedFactors=zeros(n_f,f_o);
    for q=1:f_o
        AlignedFactors(:,q)=circshift(F(:,q),-dini{w}(q));
    end
    
    Ft{w}=AlignedFactors;
    if f_o*p>1e4
        A{w}=sprandn(f_o,p,.1);    
    else
        A{w}=randn(f_o,p);
    end
    M{w}=AlignedFactors(epsilon1+1:epsilon1+n,:)*A{w};
    Noise=1/sqrt(n*p)*sqrt(norm(M{w},'fro')^2/(10^(SNR/10)))*randn(n,p);  

    MissingTP=(rand(n_f,1)>=MissingTPRatio);

    Masks{w}=sparse(n_f,p);
    X{w}=sparse(n_f,p);
    
    Masks{w}(MissingTP,:)=ones(length(find(MissingTP)),p);
    X{w}(MissingTP,:)=M{w}(MissingTP,:)+Noise(MissingTP,:);
    
    SNRest=SNRest+norm(M{w},'fro')/norm(Noise,'fro');
end
SNRest=20*log10(1/S*SNRest);
if FigOnOff
    h4=figure('Name','Generated Data');
    for i=1:S
        MissingTP=(Masks{i}(:,1)==1);
        figure(h4)
        subplot(2,ceil(S/2),i)
        plot(find(MissingTP==1),X{i}(MissingTP,:))
        title(['Subject ',num2str(i),'/',num2str(S)]);
        xlabel('Time')
        ylabel('Magnitude')
    end
end

Sigma=cell(S,1);
SigmaAv=zeros(n_f);
for i=1:S
    Sigma{i}=cov(X{i}');
    SigmaAv=SigmaAv+1/S*Sigma{i};
end

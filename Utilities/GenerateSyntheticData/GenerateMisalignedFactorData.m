function [A_o,A_no,U,F,M_no,X,M,dini,SNR,Ft]=GenerateMisalignedFactorData(p,n,n_f,S,f_o,f_no,sigma,sigmaF,Amp,VarDelays,FigOnOff,SmoothOnOff)

epsilon=n_f-n;
epsilon1=round(epsilon/2);
epsilon2=epsilon-epsilon1;


HalfIndx=floor(n_f/4);
Rest=n_f-2*HalfIndx;
MotifDictionary=zeros(n_f,4);
MotifDictionary(:,1)=[zeros(HalfIndx,1);Amp*ones(Rest,1);zeros(HalfIndx,1)]; % upregulation
MotifDictionary(:,2)=[1.5*Amp*ones(HalfIndx,1);1.5*ones(Rest,1);ones(HalfIndx,1)]; % downregulation
MotifDictionary(:,3)=[3*ones(HalfIndx-3,1);3*ones(1,1)+3*Amp*triang(1);3*ones(Rest+2,1);ones(HalfIndx,1)];% upregulation peak
MotifDictionary(:,4)=Amp*[ones(HalfIndx+2,1);ones(1,1)-triang(1);ones(Rest-3,1);ones(HalfIndx,1)];% downregulation peak
if SmoothOnOff
    for i=1:4
        MotifDictionary(:,i)=smoothSig(MotifDictionary(:,i)',6,'triang');
    end
end
    
    
%% Generate motif-gene association
Pattern=zeros(p,f_o);
LeftIndx=1:p;
for i=1:f_o
    RandIndx=randperm(length(LeftIndx));
    if i<f_o
        RandHalf=ceil((length(LeftIndx)-1)*.7*rand(1));
    else
        RandHalf=length(LeftIndx);
    end
    CurrIndx=LeftIndx(1:RandHalf);
    LeftIndx=setdiff(LeftIndx,CurrIndx);
    Pattern(CurrIndx,i)=ones(length(CurrIndx),1);
end

%% Generate Factors
F=zeros(n_f,f_o);
F_=zeros(n_f*f_o,f_o);
for i=1:f_o
    F(:,i)=abs(MotifDictionary(:,i)+ sigmaF*randn(n_f,1));
    F(:,i)=1/norm(F(:,i))*F(:,i);
    F_(((i-1)*n_f+1):(i)*n_f,i)=F(:,i);
end
   
% Non-Ordered factors
M_no=zeros(n,f_no);
if f_no>0
    M_no=MotifDictionary(epsilon1+1:epsilon1+n,3);
    M_no=1/norm(M_no)*M_no; 
end
if FigOnOff
    figure
    subplot(3,1,1)
    imagesc(Pattern)
    subplot(3,1,2)
    plot(F)
    subplot(3,1,3)
    plot(M_no)
    cols=hsv(f_o);
    figure('Name','Factors')
    for i=1:f_o
        plot(F(:,i),'Linewidth',2,'Color',cols(i,:));hold on;
        xlabel('Time')
        legstr{i}=strcat('Factor: ',num2str(i));
    end
    legend(legstr)
end

%% Generate subject's data
SNR=0;
for w=1:S
    dini{w}=floor(rand(f_o,1)*(sqrt(12*VarDelays+1)))+1;
    dini{w}=sort(abs(dini{w}),'ascend');
    for i=1:f_o
        if dini{w}(i)==0
            dini{w}(i)=1;
        end
    end
    U{w}=[];
    t=0:n-1;
    AlignedFactors=zeros(n_f,f_o);
    %Showd=dini{w}
    for q=1:f_o
        [D]=CreateTranslationMatrix(dini{w}(q),n_f);
        AlignedFactors(:,q)=D*F(:,q);
        U{w}=[U{w} D];
    end
    Ft{w}=AlignedFactors;
    A_o{w}=abs(Pattern');%ones(f_o,p);
    p1=size(A_o{w}(abs(Pattern')==1),1);
    p2=size(A_o{w}(abs(Pattern')==1),2);
    A_o{w}(abs(Pattern')==1)=abs(A_o{w}(abs(Pattern')==1)+1/5*randn(p1,p2));
    A_no{w}=abs(ones(f_no,p));%ones(f_no,p);
    % Ordered factors
%     figure
%     sub(2,1,1)
%     plot(AlignedFactors)
%     subplot(2,1,2)
%     imagesc(U{w})
    if f_no>0
        Data=AlignedFactors(epsilon1+1:epsilon1+n,:)*A_o{w}+M_no*A_no{w};
        %Data=1/max(max(Data))*Data;
        Noise=sqrt(sigma)*norm(A_o{w})*randn(n,p);
        X{w}=Data+Noise; 
        M{w}=Data;
        SNR=SNR+norm(Data,'fro')/norm(Noise,'fro');
    else
        Data=AlignedFactors(epsilon1+1:epsilon1+n,:)*A_o{w};
        %Data=1/max(max(Data))*Data;
        Noise=sqrt(sigma)*norm(A_o{w})*randn(n,p);
        Data=Data-repmat(mean(Data,2),1,p);
        X{w}=Data+Noise; 
        M{w}=Data;
        SNR=SNR+norm(Data,'fro')/norm(Noise,'fro');
    end
  % X{w}=1/max(max(X{w}))*X{w};
end
SNR=20*log10(1/S*SNR);
if FigOnOff
    figure('Name','Observed Data')
    for i=1:S
        subplot(2,ceil(S/2),i)
        plot(X{i})
    end
end

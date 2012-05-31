function [A_o,A_no,U,F,M_no,X,M,dini,SNR,Ft]=GenerateMisalignedPeriodicFactorData_SignalPower(p,n,n_f,S,f_o,f_no,Theta,sigmaF,Amp,VarDelays,FigOnOff,SmoothOnOff)

epsilon=n_f-n;
epsilon1=round(epsilon/2);
epsilon2=epsilon-epsilon1;
DFT=1/sqrt(n)*dftmtx(n);
DFTinv=conj(DFT);

% HalfIndx=round(n_f*.4);%n_f-3;
% Rest=n_f-2*HalfIndx;
% MotifDictionary=zeros(n_f,4);
% MotifDictionary(:,1)=[ones(HalfIndx,1);Amp*ones(Rest,1);ones(HalfIndx,1);]; % upregulation
% %MotifDictionary(:,2)=[1.5*Amp*ones(HalfIndx,1);1.5*ones(Rest,1);ones(HalfIndx,1);]; % downregulation
% %MotifDictionary(:,3)=[3*ones(HalfIndx-3,1);3*ones(1,1)+3*Amp*triang(1);3*ones(Rest+2,1);ones(HalfIndx,1);];% upregulation peak
% %MotifDictionary(:,4)=Amp*[ones(HalfIndx+2,1);ones(1,1)-triang(1);ones(Rest-3,1);ones(HalfIndx,1);];% downregulation peak
% if SmoothOnOff
%     for i=1:4
%         MotifDictionary(:,i)=smoothSig(MotifDictionary(:,i)',4,'triang');
%     end
% end
%     
    
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
    clear h;
    h=sort(abs(randn(round(n_f/2)+1,1)));
    h=zeros(round(n_f/2)+1,1);
    h(round(.5*end):end)=1;
    h=[h;flipud(h(2:end-1))];
    F(:,i)=h;
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
         D=diag(-exp(1i*2*pi/n*[0:n-1]*dini{w}(q)));
        % Fadj(:,j)=D*DFT*F(:,j);
        %[D]=CreateTranslationMatrix(dini{w}(q),n_f);
        AlignedFactors(:,q)=D*DFT*F(:,q);
        U{w}=[U{w} D];
    end
    Ft{w}=AlignedFactors;
%     A_o{w}=abs(Pattern');%ones(f_o,p);
%     p1=size(A_o{w}(abs(Pattern')==1),1);
%     p2=size(A_o{w}(abs(Pattern')==1),2);
%    A_o{w}(abs(Pattern')==1)=abs(A_o{w}(abs(Pattern')==1)+1/5*randn(p1,p2));
%    A_no{w}=abs(ones(f_no,p));%ones(f_no,p);
    A_o{w}=randn(f_o,p);
 %A_o{w}=1/norm(A_o{w})*A_o{w};
    A_no{w}=randn(f_no,p);
 % A_no{w}=1/norm(A_no{w})*A_no{w};
    % Ordered factors
%     figure
%     sub(2,1,1)
%     plot(AlignedFactors)
%     subplot(2,1,2)
%     imagesc(U{w})
    if f_no>0
        Data=AlignedFactors(epsilon1+1:epsilon1+n,:)*A_o{w}+M_no*A_no{w};
        %Data=1/max(max(Data))*Data;
        Noise=1/sqrt(Theta)*1/sqrt(n)*crandn(n,p);
        X{w}=Data+Noise; 
        M{w}=Data;
        SNR=SNR+norm(Data,'fro')/norm(Noise,'fro');
    else
        Data=AlignedFactors(epsilon1+1:epsilon1+n,:)*A_o{w};
        %Data=1/max(max(Data))*Data;
        Noise=1/sqrt(Theta)*1/sqrt(n)*crandn(n,p);
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
    k=1;
    for i=1:S
        subplot(S,2,k)
        plot(abs(X{i}))
        title('Modulus')
        subplot(S,2,k+1)
        plot(angle(X{i}))
        title('Phase')
        k=k+2;
    end
    figure('Name','Generated Factor')
    k=1;
    for i=1:S
         D=diag(exp(-1i*2*pi/n*[0:n-1]*dini{i}(1)));
        subplot(S,2,k)
        plot(abs(DFTinv*D*DFT*F(:,1)))
        title('Modulus')
        subplot(S,2,k+1)
        plot(angle(DFTinv*D*DFT*F(:,1)))
        title('Phase')
        k=k+2;
    end
end
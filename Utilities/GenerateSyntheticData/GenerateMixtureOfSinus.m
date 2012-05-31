function [Xsampled,XsampledF,Data,MisData,NoiseLessXd,NoiseLessAligned,SamplingTimes,driftValDiscrete]=GenerateMixtureOfSinus(S,Amp,n,N,OverSampling,k,NoiseVar)

f=1;
w=2*pi*f/N;
SamplingTimes=round(1:N/(k*OverSampling*100*f):N);
driftVal=zeros(S,1);
%NoiseLessAligned=sqrt(Amp)*1/sqrt(2)*randn*(sin(w*(n))+sin(w*3*(n ))+3*sin(w*1.5*(n )+2*sin(w*2.8*(n ))));
NoiseLessAligned=sqrt(Amp)*1/sqrt(2)*randn*(sin(w*(n))+sin(w*3*(n ))+cos(w*100*(n )));
DFT=1/sqrt(length(SamplingTimes))*dftmtx(length(SamplingTimes));
DFTinv=conj(DFT);
Xd=cell(S,1);
X=cell(S,1);
NoiseLessXd=cell(S,1);
Data=[];
MisData=[];
for i=1:S
    driftVal(i)=sign(1/2-rand)*round(N/2*rand(1));
    NoiseLessXd{i}=circshift(NoiseLessAligned,driftVal(i));
    Xd{i}=(NoiseLessXd{i}+sqrt(NoiseVar)*randn(1,N))'; %+
    X{i}=(sqrt(Amp)*randn*(NoiseLessAligned)+sqrt(NoiseVar)*randn(1,N))';
    Xsampled{i}=Xd{i}(SamplingTimes);
    XsampledF{i}=DFT*Xsampled{i};
    Data=[Data X{i}(SamplingTimes)];
    MisData=[MisData Xd{i}(SamplingTimes)];
end

driftValDiscrete=round(driftVal/(N/(k*4*OverSampling*3*f)));
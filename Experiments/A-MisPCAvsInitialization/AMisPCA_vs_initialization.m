%% FMisPCA vs Initialization

addpath(genpath('../../Utilities'))
more off;
clear all;

if ~exist("OCTAVE_VERSION","builtin")
	if matlabpool('size') == 0
	    matlabpool open 4;
	end
end

Nrndm=4*5;
NpointsSNR=4;
NpointsDmax=3;
thetaVec=logspace(-1,1,NpointsSNR);
Spacing=1;
p=100;
n=50;
r=p/10;
S=5;
F=2;

Masks=cell(S,1);
for i=1:S
    Masks{i}=ones(p,1);
end 
    
dmaxVec=round(linspace(1,p-1,NpointsDmax));
SNR=round(linspace(-10,30,NpointsSNR));

[dummy,dummy1,H]=GenerateFMisPCAData(p,n,1,F,SNR(1),r,1,Spacing,[],[]);


theta=1;
disp(sprintf('MisPCA versus Initialization, theta= %d \n',theta))

for l=1:NpointsDmax
    dmax=dmaxVec(l)       
    for k=1:NpointsSNR
	fprintf('\t dmax=%d/%d, SNR=%d/%d \n',l,NpointsDmax,k,NpointsSNR)               
        % Parallelizable for (deactivated for OCTAVE compatibility) 
	 for i=1:Nrndm
            % Generate misaligned signals with random delays
            [X,SigmaAv,Fo,dini,Sigma,D]=GenerateFMisPCAData(p,n,S,F,SNR(k),r,dmax,Spacing,[],H);
            
            Options=Simulations_options(p,max(round(dmax/4),1),dmax);
            Options.k=F;
            nvec=randn(p,F);
            hini=theta*Fo+(1-theta)*nvec; 
            [Options.Hini,dummy]=eigs(hini*hini',F);
            
            [F_BFMisPCA,dhat,Lambdas]=BF_FMisPCA(Sigma,Options);

            F_AMisPCA=FMisPCA(Sigma,Masks,Options);       
            
            CorrMisPCA(k,i,l)=MinFactorsSubspaceDistance(Fo,F_AMisPCA);   
            
            CorrBFMisPCA(k,i,l)=MinFactorsSubspaceDistance(Fo,F_BFMisPCA); 
        end
    end
end

meanAMisPCA1=squeeze(mean(CorrMisPCA,2))
stdAMisPCA1=squeeze(std(CorrMisPCA,0,2));

meanBFMisPCA=squeeze(mean(CorrBFMisPCA,2))
stdBFMisPCA=squeeze(std(CorrBFMisPCA,0,2));

%% theta=1/2;
theta=1/2;
fprintf('MisPCA versus Initialization, theta= %d \n',theta)
for l=1:NpointsDmax

    dmax=dmaxVec(l)
       
    for k=1:NpointsSNR
        fprintf('\t dmax=%d/%d, SNR=%d/%d \n',l,NpointsDmax,k,NpointsSNR)
               
        % Parallelizable for (deactivated for OCTAVE compatibility) 
	 for i=1:Nrndm
            % Generate misaligned signals with random delays
            [X,SigmaAv,Fo,dini,Sigma,D]=GenerateFMisPCAData(p,n,S,F,SNR(k),r,dmax,Spacing,[],H);
            
            Options=Simulations_options(p,max(round(dmax/4),1),dmax);
            Options.k=F;
            nvec=randn(p,F);
            hini=theta*Fo+(1-theta)*nvec; 
            [Options.Hini,dummy]=eigs(hini*hini',F);
            
            F_AMisPCA=FMisPCA(Sigma,Masks,Options);       
            
            CorrMisPCA2(k,i,l)=MinFactorsSubspaceDistance(Fo,F_AMisPCA);   
            
        end
    end
end

meanAMisPCA2=squeeze(mean(CorrMisPCA2,2))
stdAMisPCA2=squeeze(std(CorrMisPCA2,0,2));

%% theta=0;
theta=0;
fprintf('MisPCA versus Initialization, theta= %d \n',theta)
for l=1:NpointsDmax

    dmax=dmaxVec(l)
       
    for k=1:NpointsSNR

        fprintf('\t dmax=%d/%d, SNR=%d/%d \n',l,NpointsDmax,k,NpointsSNR)
               
        % Parallelizable for (deactivated for OCTAVE compatibility) 
	 for i=1:Nrndm
            % Generate misaligned signals with random, delays
            [X,SigmaAv,Fo,dini,Sigma,D]=GenerateFMisPCAData(p,n,S,F,SNR(k),r,dmax,Spacing,[],H);
            
            Options=Simulations_options(p,max(round(dmax/4),1),dmax);
            Options.k=F;
            nvec=randn(p,F);
            hini=theta*Fo+(1-theta)*nvec; 
            [Options.Hini,dummy]=eigs(hini*hini',F);
            
            F_AMisPCA=FMisPCA(Sigma,Masks,Options);       
            
            CorrMisPCA3(k,i,l)=MinFactorsSubspaceDistance(Fo,F_AMisPCA);   
            
        end
    end
end

meanAMisPCA3=squeeze(mean(CorrMisPCA3,2))
stdAMisPCA3=squeeze(std(CorrMisPCA3,0,2));


%% Latex conversion

s = MeanSTDLatexTable([SNR',meanAMisPCA1],stdAMisPCA1,'%.2g')
s2 = MeanSTDLatexTable([SNR',meanAMisPCA2],stdAMisPCA2,'%.2g')
s3 = MeanSTDLatexTable([SNR',meanAMisPCA3],stdAMisPCA3,'%.2g')
s4 = MeanSTDLatexTable([SNR',meanBFMisPCA],stdBFMisPCA,'%.2g')

fprintf('MisPCA versus Initialization: saving results... \n',theta)
save(['Results/MisPCA_vs_Initialization_',date,'.mat']);


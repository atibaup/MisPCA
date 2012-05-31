function [Options]=PHD_Analysis_Sx_Options(name,SxOrASx,p,f,n,Padding,dmax)
Options.n=n; %% #timepoints
if Padding==0
    Options.n_f=Options.n+dmax;
else
    Options.n=Options.n+2*Padding;
    Options.n_f=Options.n;
end
Options.dmax=dmax;
Options.Padding=Padding;

Options.p=p; %  #genes

Options.Nrndm=1;
Options.NpointsLambda=1; % Loadongs penalty
Options.NpointsBeta=1; % Factors penalty
Options.Ncv=1;
Options.Ntest=0; %number of subjects used for CV error estimation
Options.FiguresOnOff=0;
Options.display=1;
Options.WarmIni=1; % Compute regularization grid using warm starts
Options.Parallelize=1; % Parallalelize loop in regularization grid
Options.Dtot=eye(Options.n_f)+diag(-ones(Options.n_f-1,1),-1);
Options.GLassoOn=1;
Options.SedumiOnOff=0;
Options.PositiveLoadings=1;
Options.PositiveFactors=1;
Options.L1orL2TV=2; % L1 or L2 total variation penalty on the factors
Options.RndmOnOff=0; % "1" Random inizialization (to use with BCDMatrixRecovery_v2) or "0" (to use with BCDMatrixRecovery_v3)
Options.PenalizeMno=1;
Options.A_o=[];

subjectIDs=[1:17];
if SxOrASx==1
    Options.S=9;
    Options.SxSubjectIDs=[1 5 6 7 8 10 12 13 15];
else
    Options.S=8;
    Options.SxSubjectIDs=setdiff(1:17,[1 5 6 7 8 10 12 13 15]);
end
Options.p=p; %99;%% #genes
Options.f_o=f;  % #ordered factors
Options.f_no=0;  % #non -ordered factors
Options.Niter=5;
Options.Dtot=eye(Options.n)+1/2*diag(-ones(Options.n-1,1),-1)+1/2*diag(-ones(Options.n-1,1),1);
Options.Dtot=eye(Options.n)+diag(-ones(Options.n-1,1),-1);
Options.Initialization=3; % 1: from motif 
                          % 2: random from data 
                          % 3: from initial clusters with pseudo inverse for loadings 
                          % 4: initial clusters with hard assignemnt to loadings
Options.RndmOnOff=1; % Random initialization (overrides Options.Initialization )
Options.ConstrainOrPenalty=1; % 1: constrain on energy, 2: energy penalty on factors
Options.L1orL2TV=2;
Options.TolBCD=1e-5;
Options.PenalizeMno=0;

if Options.ConstrainOrPenalty==1
    Options.delta=1e2; % 1e-3; %% Factors' energy penalty or Factors energy constraint
else
    Options.delta=1e-4;
end
Options.BetaMin=1e-7;
Options.BetaMax=1e-3;
Options.LambdaMin=0;
Options.LambdaMax=1e-5;
Options.Lambda= 1e-30;%4/(Options.p*Options.n_f);%1e-10;%1e-9; %% Loadings' penalty
Options.Beta=1e-30; %% Factors' penalty

if strcmp(cellstr(name),'PC-de-Arnau')
    Options.NpointsBeta=1; % Factors
    Options.NpointsLambda=1;  % Loadings
    Options.Nrndm=3;%; 8*2; % Nservers*integer 
    Options.Niter=10;
else
    Options.NpointsBeta=6;
    Options.NpointsLambda=1;
    Options.Nrndm=8;%; 8*2; % Nservers*integer 
end
Options.OneByOne=100;% increasing size of the active set at each iteration
Options.display=1;
Options.FiguresOnOff=0;
Options.Parallelize=0;
Options.WarmIni=0; % '0': no warm initialization from previous solutions, '1': Warm initialization '2': warm ini only for factors 
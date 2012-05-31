function Options=DMR_PHD_DefaultOptions(SxOrASx,y,n)

Options.n=n; %% #timepoints
Options.n_f=n;
subjectIDs=[1:17];%setdiff([1:17],[13,8,4,17]); 
%SxSubjectIDs=subjectIDs(find(y==0))

Options.p=300; %99;%% #genes

Options.Mode='Delay'; % 'Delay' /'Phase'
Options.delta=1; % Energy penalty on factors
Options.Niter=20;
Options.TolBCD=1e-6;
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
Options.Initialization=2; % '1' initialize from motif dictionary, '2' initialize from data 'average'
Options.delayEstimation=1; % '1' Approximation '2' exact!
Options.A_o=[];
if SxOrASx==1
        %% Sx Analysis
    Options.f_o=1;  % #ordered factors
    Options.f_no=1;  % #non -ordered factors
    Options.Lambda=1e-2; %%   Loadings' penalty
    Options.Beta=5e-2; %%Factors' penalty 
    Options.BetaMin=Options.Beta;%1e-3;
    Options.BetaMax=Options.Beta;%1e-3;
    Options.LambdaMin=Options.Lambda;
    Options.LambdaMax=Options.Lambda;
elseif SxOrASx==2
    %% Asx Analysis
    Options.f_o=2;  % #ordered factors
    Options.f_no=1;  % #non -ordered factors
    Options.Lambda=1e-1;%% Factors' penalty 
    Options.Beta=1e-2;%% Loadings' penalty
    Options.BetaMin=Options.Beta;%1e-3;
    Options.BetaMax=Options.Beta;%1e-3;
    Options.LambdaMin=Options.Lambda;
    Options.LambdaMax=Options.Lambda;
else
    %% Both
    Options.f_o=3;  % #ordered factors
    Options.f_no=1;  % #non -ordered factors
    Options.Lambda=5e-1; %% Factors' penalty 
    Options.Beta=5e-2; %% Loadings' penalty
    Options.BetaMin=Options.Beta;%1e-3;
    Options.BetaMax=Options.Beta;%1e-3;
    Options.LambdaMin=Options.Lambda;
    Options.LambdaMax=Options.Lambda;
    
end
if nargin ==2
    if SxOrASx==1
        Options.SxSubjectIDs=subjectIDs(find(y==0));
    elseif SxOrASx==2
        Options.SxSubjectIDs=subjectIDs(find(y==1));
    else
        Options.SxSubjectIDs=subjectIDs;
    end
    Options.S=length(Options.SxSubjectIDs); %% #subjects
else
    Options.SxSubjectIDs=subjectIDs;
end
Options.Virus='Z'; % Virus under scope
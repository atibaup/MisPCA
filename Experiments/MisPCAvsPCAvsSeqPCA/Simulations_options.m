function [Options]=Simulations_options(p,Spacing,dmax)
Options.k=1;
Options.NpointsBeta=1;
Options.Nrndm=4;
Options.GridSize=Spacing;
Options.FiguresOnOff=0;
Options.dmax=dmax;
Options.p=p;
Options.OP_OnOff=0; % enforce order-preserving or not
Options.FixedDelaysOnOff=0;
Options.beta=0;
Options.Hini=[];
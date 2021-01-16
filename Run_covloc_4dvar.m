clear all;
%Evaluation of covariance localization for the proposed flow-dependent
%4D-Var method

for(b=1:4)
    params=struct;
    params.mJ=b;
    params.kass=60*9;
    params.nbblocksass=60*9;
    params.interpolin=1;
    params.neighborstrunc=3;
    params.pcgsteps=100;
    params.hybridc=0;
    params.prioruv=0.001;
    params.priorh=0.001;
    ret=shallowwater_4dvarnew_covloc('synthetic_21_10_days.mat',params);
    save(strcat("results_covloc_4dvar",num2str(b),".mat"));
end

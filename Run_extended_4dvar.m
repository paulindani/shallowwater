clear all;

%Evaluation of finding a better minimizer by gradually extending the 4D-Var
%assimilation window, as proposed by Pires et al.

arr_T=[3,6,9,12,18];

for it=1:5
    params=struct;
    params.kass=60*arr_T(it);
    params.nbblocksass=60*arr_T(it);
    params.interpolin=1;
    params.neighborstrunc=3;
    params.pcgsteps=100;
    params.mJ=0;
    params.prioruv=0.001;
    params.priorh=0.001;
    params.mxextend=2;
    ret=shallowwater_4dvarnew_extend('synthetic_21_10_days.mat',params);
    save(strcat("results_extended_4dvar_",num2str(it),".mat"));
end



%Comparison of the proposed flow-dependent 4D-Var method and Hybrid 4D-Var
%Syntethic data, 1 day experiment
%Proposed flow-dependent 4D-Var
clear all;

for(b=0:5)
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
    ret=shallowwater_4dvarnew('synthetic_21_10_days.mat',params);
    save(strcat("results_synthetic_10_days_4dvar",num2str(b),".mat"));
end
%relative errors are accessible in the ret.reldiffuv and ret.reldiffh
%components

clear all;

%Hybrid 4D-Var with various parameter choices


arr_loc=[2,3,4]
arr_inf=[5e-4,1e-3,2e-3];
arr_hyb=[0,0.1,0.3,0.5,0.7,0.9]

for(i1=1:3)
    for(i2=1:3)
        for(i3=1:6)
            params=struct;
            
            params.loc=arr_loc(i1);
            params.inflation=arr_inf(i2);
            params.hybridc=arr_hyb(i3);

            params.kass=60*9;
            params.nbblocksass=60*9;
            params.interpolin=1;
            params.neighborstrunc=3;
            params.pcgsteps=100;
            params.Nparticles=200;
            params.mJ=0;
            params.prioruv=0.001;
            params.priorh=0.001;
            
            ret=shallowwater_hybrid4dvar_climat('synthetic_21_10_days.mat',params);
            save(strcat("results_synthetic_10_days_hybrid",num2str(i1),"_",num2str(i2),"_",num2str(i3), ".mat"));
        end
    end
end


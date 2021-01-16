clear all;
%Comparison of the proposed flow-dependent 4D-Var method and Hybrid 4D-Var
%Tsunami data, 30 minutes experiment

%Proposed flow-dependent 4D-Var

for(b=0:2)
    params=struct;   
    params.mJ=b;
    params.kass=30;
    params.nbblocksass=30;
    params.interpolin=1;
    params.neighborstrunc=2;
    params.pcgsteps=200;
    params.prioruv=0.001;
    params.priorh=0.001;
    ret=shallowwater_4dvarnew('tsunami_336_30_mins.mat',params);
    save(strcat("results_tsunami_4dvar_",num2str(b),".mat"));
end

clear all;


%Hybrid 4D-Var method 

arr_loc=[2,3,4]
arr_inf=[0.01,0.02,0.03];
arr_hyb=[0,0.1,0.3,0.5,0.7,0.9]

for(i1=1:3)
    for(i2=1:3)
        for(i3=1:6)
            params=struct;
            
            params.loc=arr_loc(i1);
            params.inflation=arr_inf(i2);
            params.hybridc=arr_hyb(i3);

            params.kass=30;
            params.nbblocksass=30;
            params.interpolin=1;
            params.neighborstrunc=2;
            params.pcgsteps=200;
            params.Nparticles=100;
            params.mJ=0;
            params.prioruv=0.001;
            params.priorh=0.001;
            
            ret=shallowwater_hybrid4dvar('tsunami_336_30_mins.mat',params);
            save(strcat("results_tsunami_hybrid",num2str(i1),"_",num2str(i2),"_",num2str(i3), ".mat"));
        end
    end
end



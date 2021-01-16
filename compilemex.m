global g;
global Delta;
global f;
global cb;
global nu;
%global Blist1;
%global Blist2;
%global Blist4;
%Blist4=0;Blist1=0;Blist2=0;
g=9.81;Delta=10000;f=1e-4;cb=1e-5;nu=1e-3;

codegen dercalc.m -o dercalc_mex
codegen Dvecgen.m -o Dvecgen_mex
codegen genAmxlist.m -o genAmxlist_mex
codegen genBlist4.m -o genBlist4_mex
codegen gibbssolve1.m -o gibbssolve1_mex
codegen Jvvectomx.m -o Jvvectomx_mex
codegen localizedcov.m -o localizedcov_mex
codegen localizedcov2.m -o localizedcov2_mex
codegen neighbors_listmaker.m -o neighbors_listmaker_mex
codegen onestep.m -o onestep_mex
%codegen neighbors_listmaker0.m -o neighbors_listmaker0_mex
%codegen neighbors_listmaker2.m -o neighbors_listmaker2_mex
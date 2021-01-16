function u0it=pushforward(dermax, nbblocks, u0, Tmax,d,H,Hip1,Him1,Hjp1,Hjm1)%,BArray,FArray,ddd)
global BinomialArray;
global FactorialArray;
d23=3*d^2;
uT=zeros(2,d23);
u0it=u0;
for it=1:nbblocks
uT=onestep_mex(u0it, [0,Tmax/nbblocks], 1,dermax,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray,FactorialArray);
u0it=uT(2,1:d23);
end
end


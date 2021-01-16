function gsmval=gsmwithprioren4dvar(dermax,nbblocks,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)

global priormean;
global C;

gsmvalprior=(v-priormean)*pcg(C,(v-priormean)',1e-10,500)/2;
gsmval=gsm(dermax,nbblocks,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)+gsmvalprior;
end
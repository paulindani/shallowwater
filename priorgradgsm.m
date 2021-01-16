function [gsmval,grad]=priorgradgsm(v,priorprecmx,priormean)

gsmval=((v-priormean)*priorprecmx*(v-priormean)')/2;
 
grad=(v-priormean)*priorprecmx;
 

end

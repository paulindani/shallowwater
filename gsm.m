% function gsmval=gsm(imax,v,Y,T,k)
% global d;
% global dobs;
% global FactorialArray;
% 
% d23=3*d^2;
% 
% gsmval=0;
% Deltaobs=zeros(k+1,dobs);
% 
% vT = zeros(k+1,d23);
% vT(1, 1:d23) = v(1:d23);
% 
% vder=dercalc(imax, v);
% 
% for it = 1:k
%   vT(it+1,1:d23)=v(1:d23);
%   for j = 1:imax
%      vT(it+1,1:d23)=vT(it+1,1:d23)+vder(j + 1,1:d23)*(T(it + 1))^j/FactorialArray(j + 1);
%   end
% end
% 
% PhiT=Hfuncmx(vT,k);
% Deltaobs=Y-PhiT;
% 
% gsmval=sum(sum(Deltaobs.^2))/2;
% 
% 
% end

function gsmval=gsm(dermax,nbblocks,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)
d23=d^2*3;
vTarr=runshallowwater(dermax, nbblocks, v, T, k,d,H,Hip1,Him1,Hjp1,Hjm1);
vTarrobs=Hfuncmx(vTarr);

gsmval=sum(sum((vTarrobs(2:k+1,1:d23)-Y(2:k+1,1:d23)).^2))/2;
end






 
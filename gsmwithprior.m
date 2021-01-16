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

function gsmval=gsmwithprior(dermax,nbblocks,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)
%global d;
%global priorprecmx;
global priormean;

%d2=d^2;
%d23=d2*3;

%     u2=reshape(v(1:d2),[d,d]);
%     v2=reshape(v(d2+1:2*d2),[d,d]);
%     h2=reshape(v(2*d2+1:d23),[d,d]);
    
%verr=norm(v-priormean)
gsmvalprior=(v-priormean)*HessJcallprior((v-priormean)')/2;
% if(gsmvalprior<0)gsmvalprior
% end
gsmval=gsm(dermax,nbblocks,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)+gsmvalprior;%((v-priormean)*priorprecmx*(v-priormean)')/2;
end






 
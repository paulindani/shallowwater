
function gsmval=gsmJ(v,vT,Y)
%global priorprecmx;
global priormean;

vTobs=Hfuncmx(vT);
gsmval=sum(sum((vTobs(2:end,:)-Y(2:end,:)).^2))/2+(v-priormean)*HessJcallprior((v-priormean)')/2;
%+((v-priormean)*priorprecmx*(v-priormean)')/2;

%gsmval=sum(sum((vTobs(2:end,:)-Y(2:end,:)).^2))/2;%+((v-priormean)*priorprecmx*(v-priormean)')/2;
end






 
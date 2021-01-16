
function gsmval=gsmJen4dvar(v,vT,Y)
global priorprecmx;
global priormean;
global C;
global d;
%hiba=norm(full(priorprecmx*C-speye(d^2*3)))

vTobs=Hfuncmx(vT);
gsmval=sum(sum((vTobs(2:end,:)-Y(2:end,:)).^2))/2+(v-priormean)*pcg(C,(v-priormean)',1e-10,500)/2;
%+((v-priormean)*priorprecmx*(v-priormean)')/2;
%+((v-priormean)*priorprecmx*(v-priormean)')/2;

%gsmval=sum(sum((vTobs(2:end,:)-Y(2:end,:)).^2))/2;%+((v-priormean)*priorprecmx*(v-priormean)')/2;
end






 
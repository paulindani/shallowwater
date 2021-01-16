
function grad=gradJen4dvar(v,vT,Y,Jarr,k,d)
global priorprecmx;
global priormean;
global interpol;
global C;

d23=d^2*3;
vTobs=Hfuncmx(vT);
diffY=Y-vTobs;

grad=diffY(k+1,1:d23)';
%ngradJ0=norm(grad)


for(l=k:-1:1)
%grad=diffY(l,1:d23)'+Jarr{l}'*grad;
if(l==k|| l==1|| l==2|| mod(l,interpol)==0) grad=diffY(l,1:d23)'*(l>1)+Jarr{l}'*grad;
else
    lless=max(1,l-mod(l,interpol));
    lmore=min(k,l+interpol-mod(l,interpol));
    lambda=(l-lless)/(lmore-lless);
    grad=diffY(l,1:d23)'+Jarr{lless}'*grad*(1-lambda)+Jarr{lmore}'*grad*lambda;
end
end

%grad=grad+priorprecmx*(v-priormean)';%-diffY(1,1:d23)';
grad=grad+pcg(C,(v-priormean)',1e-6,100);
%not needed
%grad=grad-diffY(1,1:d23)';
grad=-grad;
%grad=-Jarr{1}'*(diffY(2,:))';
end






 
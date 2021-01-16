function priorprecen4dvar(d)
%global d;
global Delta;
global prioruv;
global priorh;
global prior1;
global priorprecmx;
global C;
global priormean;
%global SSOARcov;
d2=d^2;d23=3*d2;


C=spalloc(d23,d23,d23);
C(1:2*d2,1:2*d2)=speye(2*d2)*(1/prioruv);
C(2*d2+1:d23,2*d2+1:d23)=speye(d2)*(1/priorh);

priorprecmx=spalloc(d23,d23,d23);
priorprecmx(1:2*d2,1:2*d2)=speye(2*d2)*prioruv;
priorprecmx(2*d2+1:d23,2*d2+1:d23)=speye(d2)*priorh;


priormean=zeros(1,d23);


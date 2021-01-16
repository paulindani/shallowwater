function Hessw=HessJ(w,Jarr,k,d)
global priorprecmx;
global priormean;
global spobsuvprec;

d23=d^2*3;
Hessw=zeros(d23,1);

warr=zeros(d23,k+1);
warr(1:d23,1)=w;

for(l=1:k)
warr(1:d23,l+1)=(Jarr{l}')*warr(1:d23,l);
end

for(l=k:-1:1)
    Hessw=Jarr{l}*spobsuvprec*warr(1:d23,l+1)+spobsuvprec*warr(1:d23,l);
end

Hessw=Hessw+priorprecmx*w;

% vTobs=Hfuncmx(vT);
% diffY=Y-vTobs;
% 
% grad=diffY(k+1,1:d23)';
% 
% for(l=k:-1:1)
% grad=diffY(l,1:d23)'+Jarr{l}*grad;
% end
% grad=grad+priorprecmx*(v-priormean)';

end






 
function Hessw=HessJcall(w)
global kcall;
global dcall;
global Jarrcall;

global priorprecmx;
global priormean;
global spobsuvprec;
global interpol;
d=dcall;k=kcall;
d23=d^2*3;

warr=zeros(d23,k);
warr(1:d23,1)=(Jarrcall{1})*w;

for(l=2:k)

if(l==k|| l==1|| l==2|| mod(l,interpol)==0) warr(1:d23,l)=(Jarrcall{l})*warr(1:d23,l-1);
else
    lless=max(1,l-mod(l,interpol));
    lmore=min(k,l+interpol-mod(l,interpol));
    lambda=(l-lless)/(lmore-lless);
    warr(1:d23,l)=(Jarrcall{lless})*warr(1:d23,l-1)*(1-lambda)+(Jarrcall{lmore})*warr(1:d23,l-1)*lambda;
end


end

Hessw=spobsuvprec*warr(1:d23,k);
for(l=k:-1:2)
    if(l==k|| l==1|| l==2|| mod(l,interpol)==0) Hessw=Jarrcall{l}'*Hessw+spobsuvprec*warr(1:d23,l-1);
    else
    lless=max(1,l-mod(l,interpol));
    lmore=min(k,l+interpol-mod(l,interpol));
    lambda=(l-lless)/(lmore-lless);        
    Hessw=Jarrcall{lless}'*Hessw*(1-lambda)+Jarrcall{lmore}'*Hessw*lambda+spobsuvprec*warr(1:d23,l-1);
    end
end
Hessw=Jarrcall{1}'*Hessw;

Hessw=Hessw+priorprecmx*w;

% d23=dcall^2*3;
% 
% warr=zeros(d23,kcall+1);
% warr(1:d23,1)=w;
% 
% for(l=1:kcall)
% warr(1:d23,l+1)=(Jarrcall{l})*warr(1:d23,l);
% end
% 
% Hessw=spobsuvprec*warr(1:d23,kcall+1);
% for(l=kcall:-1:1)
%     Hessw=Jarrcall{l}'*Hessw+spobsuvprec*warr(1:d23,l);
% end
% 
% Hessw=Hessw-spobsuvprec*w+priorprecmx*w;
end






 
function Hessw=HessJnoprior(w,Jarr,k,d)
%global Jarrcall;
global spobsuvprec;
global interpol;

d23=d^2*3;

warr=zeros(d23,k);
warr(1:d23,1)=(Jarr{1})*w;

for(l=2:k)
%warr(1:d23,l)=(Jarr{l})*warr(1:d23,l-1);

if(l==k|| l==1|| l==2|| mod(l,interpol)==0) warr(1:d23,l)=(Jarr{l})*warr(1:d23,l-1);
else
    lless=max(1,l-mod(l,interpol));
    lmore=min(k,l+interpol-mod(l,interpol));
    lambda=(l-lless)/(lmore-lless);
    warr(1:d23,l)=(Jarr{lless})*warr(1:d23,l-1)*(1-lambda)+(Jarr{lmore})*warr(1:d23,l-1)*lambda;
end


end

Hessw=spobsuvprec*warr(1:d23,k);
for(l=k:-1:2)
    %Hessw=Jarr{l}'*Hessw+spobsuvprec*warr(1:d23,l-1);
    if(l==k|| l==1|| l==2|| mod(l,interpol)==0) Hessw=Jarr{l}'*Hessw+spobsuvprec*warr(1:d23,l-1);
    else
    lless=max(1,l-mod(l,interpol));
    lmore=min(k,l+interpol-mod(l,interpol));
    lambda=(l-lless)/(lmore-lless);        
    Hessw=Jarr{lless}'*Hessw*(1-lambda)+Jarr{lmore}'*Hessw*lambda+spobsuvprec*warr(1:d23,l-1);
    end
end
Hessw=Jarr{1}'*Hessw;

end

% d23=dcall^2*3;
% 
% warr=zeros(d23,kcall+1);
% warr(1:d23,1)=w;
% 
% for(l=1:kcall)
% warr(1:d23,l+1)=(Jarr{l})*warr(1:d23,l);
% end
% 
% Hessw=spobsuvprec*warr(1:d23,kcall+1);
% for(l=kcall:-1:1)
%     Hessw=Jarr{l}'*Hessw+spobsuvprec*warr(1:d23,l);
% end
% 
% Hessw=Hessw-spobsuvprec*w+priorprecmx*w;






 
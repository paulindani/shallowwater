function Hessw=HessJcallen4dvar(w)
global kcall;
global dcall;
global Jarrcall;
global C;

global priorprecmx;
global priormean;
global spobsuvprec;

% Hessw=w+C*(HessJcall(w)-priorprecmx*w);
% return;

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
% Hessw=Hessw-spobsuvprec*w;%+priorprecmx*w;
% %Hessw=Hessw-spobsuvprec*w+priorprecmx*w;
% %computing (I+C*H)*w, as we want to write (C^-1+H)^-1 w=(I+CH)^-1 * Cw
% Hessw=w+C*Hessw;

global kcall;
global dcall;
global Jarrcall;

global C;
%global priormean;
global spobsuvprec;
global interpol;
d=dcall;k=kcall;
d23=d^2*3;

warr=zeros(d23,k);
warr(1:d23,1)=(Jarrcall{1})*w;

for(l=2:k)
%warr(1:d23,l)=(Jarrcall{l})*warr(1:d23,l-1);

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

Hessw=w+C*Hessw;

end




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
% Hessw=Hessw-spobsuvprec*w;

%computing (I+C*H)*w, as we want to write (C^-1+H)^-1 w=(I+CH)^-1 * Cw
%Hessw=w+C*Hessw;

%+pcg(C,w,1e-6,100);
%Hessw=Hessw+priorprecmx*w;

%not needed


%Hessw=Hessw;
% vTobs=Hfuncmx(vT);
% diffY=Y-vTobs;
% 
% grad=diffY(k+1,1:d23)';
% 
% for(l=k:-1:1)
% grad=diffY(l,1:d23)'+Jarrcall{l}*grad;
% end
% grad=grad+priorprecmx*(v-priormean)';

%end






 
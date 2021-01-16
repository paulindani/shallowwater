function Hessw=HessJcallnew(w)
global kcall;
global dcall;
global Jarrcall;

% global priorprecmx;
% global priormean;
% global spobsuvprec;


Hessw=HessJcallprior(w)+HessJnoprior(w,Jarrcall,kcall,dcall);

end






 
function [priormeannew,priorprecmxnew]=priorupdateinitialobs(priorprecmx,priormean,Y0)
global spobsuvprec;

priorprecmxnew=priorprecmx+spobsuvprec;

%priormeannew=(cholmod2(priorprecmxnew,priorprecmx*priormean'+spobsuvprec*Y0',-5))';
priormeannew=(pcg(priorprecmxnew,priorprecmx*priormean'+spobsuvprec*Y0',1e-6,200))';
end

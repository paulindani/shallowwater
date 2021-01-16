function SV=Sact(V)
global sigmaz;
global d;
global spobsuvprec;
%S=spobsuvprec*C*spobsuvprec+sigmaz^2*speye(d23)

SV=sigmaz^2*V+spobsuvprec*(priorcovmxact(spobsuvprec*V,d))';
end
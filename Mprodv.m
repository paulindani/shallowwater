function Mv=Mprodv(M,v,d,amax)
global indmx;
indmax=3*(1+amax*(amax+1)*2);
d2=d^2;d23=d2*3;

Mv=sum(M(1:d23,1:indmax).*v(indmx(1:d23,1:indmax)),2);

end
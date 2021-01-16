function data_writer_synthetic(d,Deltain,Tmax,k,nbblocks,sigmaz,obsfrequv,obsfreqh,obsscenario,filename)
global BinomialArray;
global FactorialArray;
global dobs;
global Delta;
global f;
global g;
global h;
global H;global Hip1; global Him1; global Hjp1; global Hjm1;
global nu;
global cb;
global obsuvprec;
global spobsuvprec;
global imaxsim;
global prioruv;
global priorh;
global prior1;
global priorprecmx;

% k=360;
% Tmax=3600;
% nbblocks=12;

Tarr=linspace(0,Tmax,k+1);

d2=d^2;d23=d2*3;

obsuvprec=zeros(1,d23);
if(obsscenario==2)

if(obsfrequv>0)
for(i=1:obsfrequv:d)
    for(j=1:obsfrequv:d)
        ij=(i-1)*d+j;
        obsuvprec(ij)=1;
        obsuvprec(d2+ij)=1;
    end
end
end

if(obsfreqh>0)
for(i=1:obsfreqh:d)
    for(j=1:obsfreqh:d)
        ij=(i-1)*d+j;
        obsuvprec(2*d2+ij)=1;
    end
end
end
end

if(obsscenario==1)
        obsuvprec(2*d2+1)=1;obsuvprec(2*d2+2)=1;
        obsuvprec(2*d2+d+1)=1;obsuvprec(2*d2+d+2)=1;
        obsuvprec(1:2*d2)=1;
end

spobsuvprec=spdiags((obsuvprec').^2,0,d23,d23);

BinomialArrayFunc(50);
FactorialArrayFunc(50);



g=9.81;
f=1e-4;
nu=1e-3;
cb=1e-5;
Delta=Deltain;
imaxsim=10;


u0=zeros(1,d23);
H=zeros(1,d2);

for i = 1:d
    for j=1:d
        ij=(i-1)*d+j;
%         u0(ij)=1+1*sin((i+j)/d*pi*2);
%         u0(d2+ij)=1-1*cos((i-j)/d*pi*2);
%         u0(2*d2+ij)=5*sin(2*i*pi/d)*cos(2*j*pi/d);
        u0(ij)=0.5+0.5*sin((i+j)/d*pi*2);
        u0(d2+ij)=0.5-0.5*cos((i-j)/d*pi*2);
        u0(2*d2+ij)=2*sin(2*i*pi/d)*cos(2*j*pi/d);
        H(ij)=100+100*(1+sin(i/d*pi*2)/2)*(1+sin(j/d*pi*2)/2);
    end
end


Hip1=shiftip1(H,d);
Him1=shiftim1(H,d);
Hjp1=shiftjp1(H,d);
Hjm1=shiftjm1(H,d);

uTarr=runshallowwater(imaxsim,nbblocks, u0, Tarr, k,d,H,Hip1,Him1,Hjp1,Hjm1);

uTobs=Hfuncmx(uTarr);

dobs=d23;

Y=uTobs+sigmaz*randn(k+1,d23).*(ones(k+1,1)*obsuvprec);

save(filename);

end
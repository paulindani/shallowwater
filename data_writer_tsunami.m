function data_writer_tsunami(Tmax,k,nbblocks,sigmaz,obsfrequv,obsfreqh,input_filename, output_filename)
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


%%%%%%%Initial definitions & choice of initial condition

%load('tsunami_data336.mat')
load(input_filename)
d=336;
%load('tsunami_data.mat')
%d=84;

Umq=0*Umq.*(H>0);
Vmq=0*Vmq.*(H>0);
h=h.*(H>0);

g=9.81;
f=1e-4;
nu=1e-3;
cb=1e-5;

imaxsim=10;

delta=1e-2;
N=3;

Tarr=linspace(0,Tmax,k+1);


d2=d^2;
d23=3*d2;

%H=1000*ones(1,d2);

H=H.*(H>=0);
H=reshape(H,[1,d2]);
Hip1=shiftip1(H,d);
Him1=shiftim1(H,d);
Hjp1=shiftjp1(H,d);
Hjm1=shiftjm1(H,d);




obsuvprec=zeros(1,d23);

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


spobsuvprec=spdiags((obsuvprec').^2,0,d23,d23);

BinomialArrayFunc(50);
FactorialArrayFunc(50);






u0=[reshape(Umq,[1,d2]),reshape(Vmq,[1,d2]),reshape(h,[1,d2])];

%Moving forward by 10 minutes from the original initial state with zero velocities
%to create the initial state for testing data assimilation algorithms

u0=pushforward(20, 20, u0, 600,d,H,Hip1,Him1,Hjp1,Hjm1);


    uTarr=runshallowwater(imaxsim,nbblocks, u0, Tarr, k,d,H,Hip1,Him1,Hjp1,Hjm1);

    uTobs=Hfuncmx(uTarr);

    dobs=d23;

    Y=uTobs+sigmaz*randn(k+1,d23).*(ones(k+1,1)*obsuvprec);

save(output_filename);

end
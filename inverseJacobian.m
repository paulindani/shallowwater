function invJ=inverseJacobian(imax,nbblocks,N,delta,v,Tass,d,H,Hip1,Him1,Hjp1,Hjm1)

% global FactorialArray;
% global BinomialArray;
% global obsuvprec;
% global priorprecmx;


d2=d^2;
d23=3*d2;



%invJ=spalloc(d23,d23,(2*N+1)^2*d23*3);

Ivecarr=zeros(1,(2*N+1)^2*d23*3);
Jvecarr=zeros(1,(2*N+1)^2*d23*3);
Vvecarr=zeros(1,(2*N+1)^2*d23*3);


r=div(d,2*N+1);
for(i=-N:N)
    for(j=-N:N)
        D=Dvecgen_mex(d,N,r,r,i,j);
        for(uvh=0:2)
        vp=v;vm=v;vpp=v;vmm=v;
        vp(uvh*d2+(1:d2))=v(uvh*d2+(1:d2))+D*delta;
        vpp(uvh*d2+(1:d2))=v(uvh*d2+(1:d2))+2*D*delta;
        vm(uvh*d2+(1:d2))=v(uvh*d2+(1:d2))-D*delta;
        vmm(uvh*d2+(1:d2))=v(uvh*d2+(1:d2))-2*D*delta;
        vpT=pushforward(imax, nbblocks, vp, -Tass,d,H,Hip1,Him1,Hjp1,Hjm1);
        vppT=pushforward(imax, nbblocks, vpp, -Tass,d,H,Hip1,Him1,Hjp1,Hjm1);
        vmT=pushforward(imax, nbblocks, vm, -Tass,d,H,Hip1,Him1,Hjp1,Hjm1);
        vmmT=pushforward(imax, nbblocks, vmm, -Tass, d,H,Hip1,Him1,Hjp1,Hjm1);
        diffT=(-vppT+8*vpT-8*vmT+vmmT)/(12*delta);
        [Ivec,Jvec,Vvec]=Jvvectomx_mex(diffT,d,0,N,r,r,i,j,uvh);
        Ivecarr(((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+1:((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+d23)=Ivec;
        Jvecarr(((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+1:((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+d23)=Jvec;
        Vvecarr(((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+1:((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+d23)=Vvec;
        %invJ=invJ+sparse(Ivec,Jvec,Vvec,d23,d23);
        end
    end
end
invJ=sparse(Ivecarr,Jvecarr,Vvecarr,d23,d23);
end

function Jpsiarrm=Jcalcnew(imax,nbblocks,N,delta,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)

d2=d^2;
d23=3*d2;
dobs=d23;


vT=runshallowwater(imax, nbblocks, v, T, k,d,H,Hip1,Him1,Hjp1,Hjm1);
PhiT=Hfuncmx(vT);
Deltaobs=PhiT(2:k+1,1:d23)-Y(2:k+1,1:d23);

gsmval=sum(sum(Deltaobs.^2))/2;


grad=zeros(1,d23);
grad2=spalloc(d23,d23,d23*(2*N+1)^2*3);

Ivecarr=zeros(k+1,(2*N+1)^2*d23*3);
Jvecarr=zeros(k+1,(2*N+1)^2*d23*3);
Vvecarr=zeros(k+1,(2*N+1)^2*d23*3);

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
        HvpT=Hfuncmx(runshallowwater(imax, nbblocks, vp, T, k,d,H,Hip1,Him1,Hjp1,Hjm1));
        HvppT=Hfuncmx(runshallowwater(imax, nbblocks, vpp, T, k,d,H,Hip1,Him1,Hjp1,Hjm1));
        HvmT=Hfuncmx(runshallowwater(imax, nbblocks, vm, T, k,d,H,Hip1,Him1,Hjp1,Hjm1));
        HvmmT=Hfuncmx(runshallowwater(imax, nbblocks, vmm, T, k,d,H,Hip1,Him1,Hjp1,Hjm1));
        HdiffT=(-HvppT+8*HvpT-8*HvmT+HvmmT)/(12*delta);
        [Ivec,Jvec,Vvec]=Jvvectomx_mex(HdiffT,d,k,N,r,r,i,j,uvh);
        Ivecarr(1:k+1,((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+1:((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+d23)=Ivec;
        Jvecarr(1:k+1,((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+1:((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+d23)=Jvec;
        Vvecarr(1:k+1,((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+1:((i+N)*(2*N+1)+(j+N))*d23+uvh*(2*N+1)^2*d23+d23)=Vvec;
         end
    end
end

for(m=1:k)
    Jpsiarrm=sparse(Ivecarr(m+1,:),Jvecarr(m+1,:),Vvecarr(m+1,:),d23,d23);
    grad=grad+((Jpsiarrm)'*Deltaobs(m,1:dobs)')';
    grad2=grad2+(Jpsiarrm)'*(Jpsiarrm);
end
%Jpsiarrm=sparse(Ivecarr(m,:),Jvecarr(m,:),Vvecarr(m,:),d23,d23);

end



 
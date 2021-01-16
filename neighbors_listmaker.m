function [Ivec,Jvec,Vvec]=neighbors_listmaker(d)
assert(isa(d, 'double'));
d2=d^2;

Ivec=zeros(15*d2,1);
Jvec=zeros(15*d2,1);
Vvec=ones(15*d2,1);
it=1;
for(i=1:d)
    for(j=1:d)
        ij=(i-1)*d+j;
        ip1j=(modd(i+1,d)-1)*d+j;
        im1j=(modd(i-1,d)-1)*d+j;
        ijp1=(i-1)*d+modd(j+1,d);
        ijm1=(i-1)*d+modd(j-1,d);
        d2ij=d2+ij;
        d2ip1j=d2+ip1j;
        d2im1j=d2+im1j;
        d2ijp1=d2+ijp1;
        d2ijm1=d2+ijm1;
        td2ij=2*d2+ij;
        td2ip1j=2*d2+ip1j;
        td2im1j=2*d2+im1j;
        td2ijp1=2*d2+ijp1;
        td2ijm1=2*d2+ijm1;
        Ivec(it)=ij;Jvec(it)=ij;it=it+1;
        Ivec(it)=ij;Jvec(it)=ip1j;it=it+1;
        Ivec(it)=ij;Jvec(it)=im1j;it=it+1;
        Ivec(it)=ij;Jvec(it)=ijp1;it=it+1;
        Ivec(it)=ij;Jvec(it)=ijm1;it=it+1;
        Ivec(it)=ij;Jvec(it)=d2ij;it=it+1;
        Ivec(it)=ij;Jvec(it)=d2ip1j;it=it+1;
        Ivec(it)=ij;Jvec(it)=d2im1j;it=it+1;
        Ivec(it)=ij;Jvec(it)=d2ijp1;it=it+1;
        Ivec(it)=ij;Jvec(it)=d2ijm1;it=it+1;
        Ivec(it)=ij;Jvec(it)=td2ij;it=it+1;
        Ivec(it)=ij;Jvec(it)=td2ip1j;it=it+1;
        Ivec(it)=ij;Jvec(it)=td2im1j;it=it+1;
        Ivec(it)=ij;Jvec(it)=td2ijp1;it=it+1;
        Ivec(it)=ij;Jvec(it)=td2ijm1;it=it+1;
    end
end

end
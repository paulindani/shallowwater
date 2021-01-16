function Amxlist=genAmxlist(d,H)

  assert(isa(d, 'double'));

  assert(isa(H, 'double'));


  assert(all(size(H)<=[1,4000000]));

%global d;
global Delta;
global f;
%global Amxlist;
global g;
global cb;
global nu;
%global H;
d2=d^2;
Amxlist=zeros(22*d2,3);
g2D=g/(2*Delta);
nuD2=nu/(Delta^2);
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
        Amxlist(it,:)=[ij,d2ij,-f];it=it+1;
        Amxlist(it,:)=[d2ij,ij,f];it=it+1;
        Amxlist(it,:)=[ij,td2ip1j,g2D];it=it+1;
        Amxlist(it,:)=[ij,td2im1j,-g2D];it=it+1;
        Amxlist(it,:)=[d2ij,td2ijp1,g2D];it=it+1;
        Amxlist(it,:)=[d2ij,td2ijm1,-g2D];it=it+1;
        Amxlist(it,:)=[ij,ij,cb+4*nuD2];it=it+1;
        Amxlist(it,:)=[d2ij,d2ij,cb+4*nuD2];it=it+1;
        
        Amxlist(it,:)=[ij,ip1j,-nuD2];it=it+1;
        Amxlist(it,:)=[ij,im1j,-nuD2];it=it+1;
        Amxlist(it,:)=[ij,ijp1,-nuD2];it=it+1;
        Amxlist(it,:)=[ij,ijm1,-nuD2];it=it+1;
        
        Amxlist(it,:)=[d2ij,d2ip1j,-nuD2];it=it+1;
        Amxlist(it,:)=[d2ij,d2im1j,-nuD2];it=it+1;
        Amxlist(it,:)=[d2ij,d2ijp1,-nuD2];it=it+1;
        Amxlist(it,:)=[d2ij,d2ijm1,-nuD2];it=it+1;
        
        Amxlist(it,:)=[td2ij,ip1j,H(ij)/(2*Delta)];it=it+1;
        Amxlist(it,:)=[td2ij,im1j,-H(ij)/(2*Delta)];it=it+1;
        Amxlist(it,:)=[td2ij,d2ijp1,H(ij)/(2*Delta)];it=it+1;
        Amxlist(it,:)=[td2ij,d2ijm1,-H(ij)/(2*Delta)];it=it+1;
        Amxlist(it,:)=[td2ij,ij,(H(ip1j)-H(im1j))/(2*Delta)];it=it+1;
        Amxlist(it,:)=[td2ij,d2ij,(H(ijp1)-H(ijm1))/(2*Delta)];it=it+1;

    end
end

end
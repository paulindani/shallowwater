function [Blist1,Blist2,Blist4]=genBlist4(d)

assert(isa(d, 'double'));
d2=d^2;


% global Blist4;
% global Blist2;
% global Blist1;
%global d;
global Delta;
%global f;
%global g;


Blist4=zeros(4*d2,4);
Blist2=zeros(12*d2,4);
Blist1=zeros(16*d2,4);


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
 
        Blist4(it,1:4)=[td2ij,td2ij,ip1j,1];it=it+1;
        Blist4(it,1:4)=[td2ij,td2ij,im1j,-1];it=it+1;
        Blist4(it,1:4)=[td2ij,td2ij,d2ijp1,1];it=it+1;
        Blist4(it,1:4)=[td2ij,td2ij,d2ijm1,-1];it=it+1;
    end
end

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
 

        Blist2(it,1:4)=[d2ij,ij,ijp1,1];it=it+1;
        Blist2(it,1:4)=[d2ij,ij,ijm1,-1];it=it+1;
        
        Blist2(it,1:4)=[ij,ij,ip1j,1];it=it+1;
        Blist2(it,1:4)=[ij,ij,im1j,-1];it=it+1;

        Blist2(it,1:4)=[ij,d2ij,d2ip1j,1];it=it+1;
        Blist2(it,1:4)=[ij,d2ij,d2im1j,-1];it=it+1;

        Blist2(it,1:4)=[d2ij,d2ij,d2ijp1,1];it=it+1;
        Blist2(it,1:4)=[d2ij,d2ij,d2ijm1,-1];it=it+1;        
        

        Blist2(it,1:4)=[ij,td2ij,td2ip1j,1];it=it+1;
        Blist2(it,1:4)=[ij,td2ij,td2im1j,-1];it=it+1;
        
        Blist2(it,1:4)=[d2ij,td2ij,td2ijp1,1];it=it+1;
        Blist2(it,1:4)=[d2ij,td2ij,td2ijm1,-1];it=it+1;
    end
end
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
        
        Blist1(it,1:4)=[ip1j,td2ij,td2ij,1];it=it+1;
        Blist1(it,1:4)=[im1j,td2ij,td2ij,-1];it=it+1;
        Blist1(it,1:4)=[d2ijp1,td2ij,td2ij,1];it=it+1;
        Blist1(it,1:4)=[d2ijm1,td2ij,td2ij,-1];it=it+1;
        
        Blist1(it,1:4)=[ijp1,ij,d2ij,1];it=it+1;
        Blist1(it,1:4)=[ijm1,ij,d2ij,-1];it=it+1;
        
        Blist1(it,1:4)=[ip1j,ij,ij,1];it=it+1;
        Blist1(it,1:4)=[im1j,ij,ij,-1];it=it+1;

        Blist1(it,1:4)=[d2ip1j,d2ij,ij,1];it=it+1;
        Blist1(it,1:4)=[d2im1j,d2ij,ij,-1];it=it+1;

        Blist1(it,1:4)=[d2ijp1,d2ij,d2ij,1];it=it+1;
        Blist1(it,1:4)=[d2ijm1,d2ij,d2ij,-1];it=it+1;        
        

        Blist1(it,1:4)=[td2ip1j,td2ij,ij,1];it=it+1;
        Blist1(it,1:4)=[td2im1j,td2ij,ij,-1];it=it+1;
        
        Blist1(it,1:4)=[td2ijp1,td2ij,d2ij,1];it=it+1;
        Blist1(it,1:4)=[td2ijm1,td2ij,d2ij,-1];it=it+1;        

    end
end
Blist4(1:end,4)=Blist4(1:end,4)/(4*Delta);
Blist2(1:end,4)=Blist2(1:end,4)/(4*Delta);
Blist1(1:end,4)=Blist1(1:end,4)/(4*Delta);
end
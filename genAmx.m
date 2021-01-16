function Amx=genAmx(d,H)
%global d;
global Delta;
global f;
%global Amx;
global g;
d23=d^2*3;d2=d^2;

%Amxlist=genAmxlist_mex;
Amxlist=genAmxlist_mex(d,H);
Amx=sparse(Amxlist(:,1),Amxlist(:,2),Amxlist(:,3),d23,d23);
%Amx=sparse(d23,d23);
% g2D=g/(2*Delta);
% for(i=1:d)
%     for(j=1:d)
%         ij=(i-1)*d+j;
%         ip1j=(modd(i+1,d)-1)*d+j;
%         im1j=(modd(i-1,d)-1)*d+j;
%         ijp1=(i-1)*d+modd(j+1,d);
%         ijm1=(i-1)*d+modd(j-1,d);
%         d2ij=d2+ij;
%         d2ip1j=d2+ip1j;
%         d2im1j=d2+im1j;
%         d2ijp1=d2+ijp1;
%         d2ijm1=d2+ijm1;
%         td2ij=2*d2+ij;
%         td2ip1j=2*d2+ip1j;
%         td2im1j=2*d2+im1j;
%         td2ijp1=2*d2+ijp1;
%         td2ijm1=2*d2+ijm1;
%         Amx(ij,d2ij)=-f;
%         Amx(d2ij,ij)=f;
%         Amx(ij,td2ip1j)=g2D;
%         Amx(ij,td2im1j)=-g2D;
%         Amx(d2ij,td2ijp1)=g2D;
%         Amx(d2ij,td2ijm1)=-g2D;
%     end
% end

end
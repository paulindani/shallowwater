function M3=neighbors(N,d)
%global d;
d2=d^2;d23=3*d2;

% Mi=[1:d23, 1:d2, 1:d2,d2+1:2*d2,d2+1:2*d2,2*d2+1:d23,2*d2+1:d23];
% Mj=[1:d23, d2+1:2*d2, 2*d2+1:d23,1:d2,2*d2+1:d23,1:d2,d2+1:2*d2];
[Ivec,Jvec,Vvec]=neighbors_listmaker_mex(d);
M=sparse(Ivec,Jvec,Vvec,d23,d23,(2*N+1)^2*3*d23);
M3=speye(d23);
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
%         M(ij,ij)=1;
%         M(ij,ip1j)=1;
%         M(ij,im1j)=1;
%         M(ij,ijp1)=1;
%         M(ij,ijm1)=1;
%         M(ij,d2ij)=1;
%         M(ij,d2ip1j)=1;
%         M(ij,d2im1j)=1;
%         M(ij,d2ijp1)=1;
%         M(ij,d2ijm1)=1;
%         M(ij,td2ij)=1;
%         M(ij,td2ip1j)=1;
%         M(ij,td2im1j)=1;
%         M(ij,td2ijp1)=1;
%         M(ij,td2ijm1)=1;
%     end
% end
M(d2+1:2*d2,1:d23)=M(1:d2,1:d23);
M(2*d2+1:d23,1:d23)=M(1:d2,1:d23);
for(it=1:N)
    M3=((M3*M)>0);
end

end
function BJvdv=BfuncJvdv(Jv,dv)
%global Blist;
global d;
%f=1e-4;Delta=1e3;g=9.81;g2D=g/(2*Delta);
d2=d^2;d23=d2*3;

%BJvdv=sparse(d23,d23);
% Bdv=sparse(d23,d23);
% I=speye(d23);
% for(i=1:d23)
% %    BJvdv(:,i)=Bfunc(Jv(:,i)',dv);
%     Bdv(:,i)=Bfunc(I(:,i)',dv);
% end
% BJvdv=Bdv*Jv;

BJvdv=Bvfunc(dv)*Jv;
%BJvdv=(BJvdv+BJvdv')/2;

end
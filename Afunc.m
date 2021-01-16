function Au=Afunc(u,d,H,Hip1,Him1,Hjp1,Hjm1)
%global d;
global Delta;
global g;
global f;
%global H;global Hip1; global Him1; global Hjp1; global Hjm1;
global cb;
global nu;

% assert(isa(d, 'double'));
% assert(isa(u, 'double'));
% assert(all(size(u)<=[1,3000000]));

g2D=g/(2*Delta);
d2=d^2;
nuD2=nu/(Delta^2);
Au=zeros(1,3*d2);
uip1=shiftip1(u(1:d2),d);
uim1=shiftim1(u(1:d2),d);
ujp1=shiftjp1(u(1:d2),d);
ujm1=shiftjm1(u(1:d2),d);
vip1=shiftip1(u(d2+1:2*d2),d);
vim1=shiftim1(u(d2+1:2*d2),d);
vjp1=shiftjp1(u(d2+1:2*d2),d);
vjm1=shiftjm1(u(d2+1:2*d2),d);

Au(1:d2)=-f*u(d2+1:2*d2)+cb*u(1:d2)-nuD2*(uip1+uim1+ujp1+ujm1-4*u(1:d2))+g2D*(shiftip1(u(2*d2+1:3*d2),d)-shiftim1(u(2*d2+1:3*d2),d));
Au(d2+1:2*d2)=f*u(1:d2)+cb*u(d2+1:2*d2)-nuD2*(vip1+vim1+vjp1+vjm1-4*u(d2+1:2*d2))+g2D*(shiftjp1(u(2*d2+1:3*d2),d)-shiftjm1(u(2*d2+1:3*d2),d));
Au(2*d2+1:3*d2)=(1/(2*Delta))*H.*(uip1-uim1+vjp1-vjm1)+(1/(2*Delta))*u(1:d2).*(Hip1-Him1)+(1/(2*Delta))*u(d2+1:2*d2).*(Hjp1-Hjm1);
%Au(2*d2+1:3*d2)=-(1/(2*Delta))*H.*(uip1-uim1+vjp1-vjm1)-(1/(2*Delta))*u(1:d2).*(Hip1-Him1)-(1/(2*Delta))*u(d2+1:2*d2).*(Hjp1-Hjm1);

end
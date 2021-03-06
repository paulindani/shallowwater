function Buv=Bfunc(ul,ur,d)
global Delta;
% assert(isa(ul, 'double'));
% assert(all(size(ul)<=[1,3000000]));
% 
% assert(isa(ur, 'double'));
% assert(all(size(ur)<=[1,3000000]));
% 

d2=d^2;
Buv=zeros(1,3*d2);


Buv(1:d2)=((-shiftjp1(ul(1:d2),d)+shiftjm1(ul(1:d2),d)).*ur(d2+1:2*d2)...
           +(-shiftjp1(ur(1:d2),d)+shiftjm1(ur(1:d2),d)).*ul(d2+1:2*d2)...
           -(shiftip1(ul(1:d2),d)-shiftim1(ul(1:d2),d)).*ur(1:d2)...
           -(shiftip1(ur(1:d2),d)-shiftim1(ur(1:d2),d)).*ul(1:d2))/(-4*Delta);
       
Buv(d2+1:2*d2)=(-(shiftip1(ul(d2+1:2*d2),d)-shiftim1(ul(d2+1:2*d2),d)).*ur(1:d2)...
           -(shiftip1(ur(d2+1:2*d2),d)-shiftim1(ur(d2+1:2*d2),d)).*ul(1:d2)...
           -(shiftjp1(ul(d2+1:2*d2),d)-shiftjm1(ul(d2+1:2*d2),d)).*ur(d2+1:2*d2)...
           -(shiftjp1(ur(d2+1:2*d2),d)-shiftjm1(ur(d2+1:2*d2),d)).*ul(d2+1:2*d2))/(-4*Delta);

Buv(2*d2+1:3*d2)=(ul(2*d2+1:3*d2).*(shiftip1(ur(1:d2),d)-shiftim1(ur(1:d2),d)+shiftjp1(ur(d2+1:2*d2),d)-shiftjm1(ur(d2+1:2*d2),d))...
                 +ur(2*d2+1:3*d2).*(shiftip1(ul(1:d2),d)-shiftim1(ul(1:d2),d)+shiftjp1(ul(d2+1:2*d2),d)-shiftjm1(ul(d2+1:2*d2),d))...
                 +ul(1:d2).*(shiftip1(ur(2*d2+1:3*d2),d)-shiftim1(ur(2*d2+1:3*d2),d))...
                 +ur(1:d2).*(shiftip1(ul(2*d2+1:3*d2),d)-shiftim1(ul(2*d2+1:3*d2),d))...
                 +ul(d2+1:2*d2).*(shiftjp1(ur(2*d2+1:3*d2),d)-shiftjm1(ur(2*d2+1:3*d2),d))...
                 +ur(d2+1:2*d2).*(shiftjp1(ul(2*d2+1:3*d2),d)-shiftjm1(ul(2*d2+1:3*d2),d)))/(4*Delta);
             
end
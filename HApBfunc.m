function HApB=HApBfunc(A, B)%, d, imax)
%   imax21=2*imax+1;
%   imax41=4*imax+1;
%   imax81=8*imax+1;
%   HApB=zeros(imax81,d);
%   for i6 = 0:div(d,6)-1
%     for iplus = 1:3
%       i1=i6*6+iplus;
%       for i2mi1 = -2*imax:2*imax
%         for i3mi1 = -2*imax:2*imax
%           HApB(imax41+i3mi1-i2mi1,modd(i1+i2mi1,d))=HApB(imax41+i3mi1-i2mi1,modd(i1+i2mi1,d))+A(imax21+i2mi1,i1)*B(imax21+i3mi1,i1);
%         end
%       end
%     end
%   end
global d;
global obsuvprec;
global spobsuvprec;
d2=d^2;
d23=d2*3;
HApB=(spobsuvprec*A)'*B;

%HApB=spdiags((obsuvprec').^2,0,d23,d23)*(A'*B);
%HApB(1:2*d2,1:d23)=HApB(1:2*d2,1:d23)*obsuvprec^2;
end


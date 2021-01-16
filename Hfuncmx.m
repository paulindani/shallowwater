function PhiT=Hfuncmx(vT)
%   PhiT=zeros(k+1,div(int64(d),int64(2)));
%   for i = 1:div(d,6)
%     PhiT(1:k+1,3*(i-1)+1)=vT(1:k+1,6*(i-1)+1);
%     PhiT(1:k+1,3*(i-1)+2)=vT(1:k+1,6*(i-1)+2);
%     PhiT(1:k+1,3*(i-1)+3)=vT(1:k+1,6*(i-1)+3);
%   end
global spobsuvprec;

% PhiT=zeros(size(vT));
% for(i=1:k+1)
% PhiT(i,1:d23)=obsuvprec.*vT(i,1:d23);
% end
PhiT=vT*spobsuvprec;
%PhiT=vT(1:(k+1),2*d2+1:d23);
end
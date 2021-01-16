function s=shiftip1(u,d)
d2=d^2;
% s=zeros(1,d2);
% s(1:(d2-d))=u((d+1):d2);
s=zeros(1,d^2);
s(1:d2-d)=u((d+1):d2);
s(d2-d+1:d2)=u(1:d);
%s=[u((d+1):d2),u(1:d)];
end
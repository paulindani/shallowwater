function s=shiftjm1(u,d)
%global d;
d2=d^2;
s=zeros(1,d2);
for(i=0:d-1)
    s((i*d+2):((i+1)*d))=u((i*d+1):((i+1)*d-1)); s(i*d+1)=u((i+1)*d);
end
end
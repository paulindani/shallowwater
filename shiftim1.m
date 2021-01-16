function s=shiftim1(u,d)
%global d;
d2=d^2;
s=zeros(1,d2);
s(1:d)=u((d2-d+1):d2);
s((d+1):d2)=u(1:(d2-d));
%s=horzcat(u((d2-d+1):d2),u(1:(d2-d)));
end
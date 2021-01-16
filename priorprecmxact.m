function PV=priorprecmxact(V,d)
%global d;
global Delta;
global prioruv;
global priorh;
global prior1;
global priorprecmx;
global priormean;
d2=d^2;d23=3*d2;

u2=reshape(V(1:d2),[d,d]);
v2=reshape(V(d2+1:2*d2),[d,d]);
h2=reshape(V(2*d2+1:d23),[d,d]);

s=0.4;
K1=(1:d)'*ones(1,d);
K2=K1';

%SSOARcov=(1+(K1.^2+K2.^2)*s^2).^(-2)*4*s;
SSOARprec=(1+(K1.^2+K2.^2)*s^2).^(2)/(4*s);

u2=idct2(SSOARprec.*dct2(u2))*prioruv;
v2=idct2(SSOARprec.*dct2(v2))*prioruv;
h2=idct2(SSOARprec.*dct2(h2))*priorh;

PV=[reshape(u2,[1,d2]),reshape(v2,[1,d2]),reshape(h2,[1,d2])];
end

% % 
% priorprecmx1=zeros(d2,d2);
% C=zeros(d2,d2);
% DperL=1;
% for(i1=1:d)
% for(j1=1:d)
% for(i2=1:d)
% for(j2=1:d)
%         di=min(abs(i1-i2),d-abs(i1-i2));
%         dj=min(abs(j1-j2),d-abs(j1-j2));
%         r=sqrt(di^2+dj^2);
%         ij1=(i1-1)*d+j1;ij2=(i2-1)*d+j2;
%         C(ij1,ij2)=(1+DperL*r)*exp(-DperL*r);
% end
% end
% end
% end
% priorprecmx1=inv(C);
% priorprecmx=zeros(d23,d23);
% priorprecmx(1:d2,1:d2)=priorprecmx1*prioruv;
% priorprecmx(d2+1:2*d2,d2+1:2*d2)=priorprecmx1*prioruv;
% priorprecmx(2*d2+1:d23,2*d2+1:d23)=priorprecmx1*priorh;





% priorprecmx=spalloc(d23,d23,d23);
% priorprecmx(1:2*d2,1:2*d2)=speye(2*d2)*prioruv;
% priorprecmx(2*d2+1:d23,2*d2+1:d23)=speye(d2)*priorh;





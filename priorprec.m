function priorprec(d)
%global d;
global Delta;
global prioruv;
global priorh;
global prior1;
global priorprecmx;
global priorcovmx;
global priormean;
%global SSOARcov;
d2=d^2;d23=3*d2;


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
% priorcovmx=sparse(d23,d23);
% priorcovmx(1:d2,1:d2)=C/prioruv;
% priorcovmx(d2+1:2*d2,d2+1:2*d2)=C/prioruv;
% priorcovmx(2*d2+1:d23,2*d2+1:d23)=C/priorh;
% priorprecmx1=inv(C);
% priorprecmx=zeros(d23,d23);
% priorprecmx(1:d2,1:d2)=priorprecmx1*prioruv;
% priorprecmx(d2+1:2*d2,d2+1:2*d2)=priorprecmx1*prioruv;
% priorprecmx(2*d2+1:d23,2*d2+1:d23)=priorprecmx1*priorh;



% 
priorprecmx=spalloc(d23,d23,d23);
priorprecmx(1:2*d2,1:2*d2)=speye(2*d2)*prioruv;
priorprecmx(2*d2+1:d23,2*d2+1:d23)=speye(d2)*priorh;
% 
% 


% priorcovmx=spalloc(d23,d23,d23);
% priorcovmx(1:2*d2,1:2*d2)=speye(2*d2)*(1/prioruv);
% priorcovmx(2*d2+1:d23,2*d2+1:d23)=speye(d2)*(1/priorh);



%priorcovmx=speye(d23);


% M=sparse(d2,d2);priorprecmx=sparse(d23,d23);
% M=speye(d2)*4;
% for(i=1:d)
%     for(j=1:d)
%         ij=(i-1)*d+j;
%         ip1j=(modd(i+1,d)-1)*d+j;
%         %im1j=(modd(i-1,d)-1)*d+j;
%         ijp1=(i-1)*d+modd(j+1,d);
%         %ijm1=(i-1)*d+modd(j-1,d);
%         M(ij,ip1j)=-1;
%         M(ip1j,ij)=-1;
%         M(ijp1,ij)=-1;
%         M(ij,ijp1)=-1;
%     end
% end
% priorprecmx(1:d2,1:d2)=M*prioruv;
% priorprecmx(d2+1:2*d2,d2+1:2*d2)=priorprecmx(1:d2,1:d2);
% priorprecmx(2*d2+1:d23,2*d2+1:d23)=M*priorh;
% if(prior1~=0)priorprecmx=priorprecmx+speye(d23)*prior1;end

priormean=zeros(1,d23);

% end
% dmax=20000/Delta;
% 
% M=sparse(d2,d2);priorprecmx=sparse(d23,d23);
% 
% for(i=1:d)
%     for(j=1:d)
%         ij=(i-1)*d+j;
%         ip1j=(modd(i+1,d)-1)*d+j;
%         im1j=(modd(i-1,d)-1)*d+j;
%         ijp1=(i-1)*d+modd(j+1,d);
%         ijm1=(i-1)*d+modd(j-1,d);
%         M(ij,ij)=1;
%         M(ij,ip1j)=1;
%         M(ij,im1j)=1;
%         M(ij,ijp1)=1;
%         M(ij,ijm1)=1;
%     end
% end
% 
% Mn=speye(d2);
% for(it=1:ceil(dmax/Delta))
%     Mn=((Mn*M)>0);
% end
% 
% clear M;
% 
% M=sparse(d2,d2);
% 
% nb=0;
% for(ij=1:d2)
%     neighborlist=find(Mn(:,ij));
%     ln=length(neighborlist);
%     for(it=1:ln)
%         ij2=neighborlist(it);
%         if(ij2>ij)
%         j1=modd(ij,d);
%         i1=(ij-j1)/d+1;
%         j2=modd(ij2,d);
%         i2=(ij2-j2)/d+1;
%         dist12=sqrt((j1-j2)^2+(i1-i2)^2);
%         c=1/dist12;
%         if(dist12>1e-10 && dist12<=(dmax+1e-10))
%         M(ij,ij)=M(ij,ij)+c;
%         M(ij2,ij2)=M(ij2,ij2)+c;
%         M(ij,ij2)=M(ij,ij2)-c;
%         M(ij2,ij)=M(ij2,ij)-c;
%         end
%         end
%     end
% end
% 
% 
% priorprecmx(1:d2,1:d2)=M*prioruv;
% priorprecmx(d2+1:2*d2,d2+1:2*d2)=priorprecmx(1:d2,1:d2);
% priorprecmx(2*d2+1:d23,2*d2+1:d23)=M*priorh;
% if(prior1~=0)priorprecmx=priorprecmx+speye(d23)*prior1;end
% end

%end
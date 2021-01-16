function Mp=trp(M,d,amax)
global inddifflist;
global invindlist;
global indmx;
indmax=size(M,2);%3*(1+amax*(amax+1)*2);
d2=d^2;d23=d2*3;

Mp=zeros(d23,indmax);
for(ind=1:indmax)
    ipdiff=inddifflist(ind,1);
    jpdiff=inddifflist(ind,2);
    for(uvh=1:3)               
    indp=invindlist(amax+1-ipdiff,amax+1-jpdiff,uvh);
    Mp(indmx((uvh-1)*d2+1:uvh*d2,ind),indp)=M((uvh-1)*d2+1:uvh*d2,ind);
    %for(i=1:d)
        %for(j=1:d)        
                 %uvhp=inddifflist(ind,3);
%                  ip=modd(i+ipdiff,d);
%                  jp=modd(j+jpdiff,d);
%                 ijuvh=(uvh-1)*d2+(i-1)*d+j;
%                  ijuvhp=(uvhp-1)*d2+(ip-1)*d+jp;             
%                 Mp(indmx(ijuvh,ind),indp)=M((uvh-1)*d2+(i-1)*d+j,ind);                
%            end
%        end
    end
end


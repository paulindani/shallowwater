function indmx=genindmx(d,amax)
global inddifflist;
indmax=3*(1+amax*(amax+1)*2);

d2=d^2;d23=d2*3;
indmx=zeros(d23,indmax);
for(i=1:d)
    for(j=1:d)
        for(uvh=1:3)
            for(ind=1:indmax)    
                 ipdiff=inddifflist(ind,1);
                 jpdiff=inddifflist(ind,2);
                 uvhp=inddifflist(ind,3);
                 ip=modd(i+ipdiff,d);
                 jp=modd(j+jpdiff,d);
                 ijuvh=(uvh-1)*d2+(i-1)*d+j;
                 ijuvhp=(uvhp-1)*d2+(ip-1)*d+jp;
                 indmx(ijuvh,ind)=ijuvhp;
        end
    end
end


end
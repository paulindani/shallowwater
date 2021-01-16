function [inddifflist,invindlist]=genindifflist(amax)
indmax=3*(1+amax*(amax+1)*2);
inddifflist=zeros(indmax,3);
invindlist=zeros(2*amax+1,2*amax+1,3);
inddifflist(1,1:3)=[0,0,1];invindlist(amax+1,amax+1,1)=1;
inddifflist(2,1:3)=[0,0,2];invindlist(amax+1,amax+1,2)=2;
inddifflist(3,1:3)=[0,0,3];invindlist(amax+1,amax+1,3)=3;
it=3;

for a=1:amax
    for ait=0:(a-1)
        for uvh=1:3
        it=it+1;
        ip=ait;jp=a-ait;
        inddifflist(it,1:3)=[ip,jp,uvh];
        invindlist(amax+1+ip,amax+1+jp,uvh)=it;
        end
    end

    for ait=0:(a-1)
        for uvh=1:3
        it=it+1;
        ip=a-ait;jp=-ait;
        inddifflist(it,1:3)=[ip,jp,uvh];
        invindlist(amax+1+ip,amax+1+jp,uvh)=it;
        end
    end
    
    for ait=0:(a-1)
        for uvh=1:3
        it=it+1;
        ip=-ait;jp=-a+ait;
        inddifflist(it,1:3)=[ip,jp,uvh];
        invindlist(amax+1+ip,amax+1+jp,uvh)=it;
        end
    end
    
    for ait=0:(a-1)
        for uvh=1:3
        it=it+1;
        ip=-a+ait;jp=ait;
        inddifflist(it,1:3)=[ip,jp,uvh];
        invindlist(amax+1+ip,amax+1+jp,uvh)=it;
        end
    end    
    
    
end
       

end
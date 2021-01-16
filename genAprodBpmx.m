function [AprodBpmxA,AprodBpmxB,AprodBpnb]=genAprodBpmx(a,b,amaxABp)
global invindlist;
%invindlist=zeros(2*amax+1,2*amax+1,3);
AprodBpmxA=zeros(2*amaxABp+1,2*amaxABp+1,3*(2*amaxABp+1)^2);
AprodBpmxB=zeros(2*amaxABp+1,2*amaxABp+1,3*(2*amaxABp+1)^2);
AprodBpnb=zeros(2*amaxABp+1,2*amaxABp+1);

for ipdiff=-amaxABp:amaxABp
    for jpdiff=(-amaxABp+abs(ipdiff)):(amaxABp-abs(ipdiff))
        
        for(idiff=max(-a,ipdiff-b):min(a,ipdiff+b))
            for(jdiff=max(-a,jpdiff-b):min(a,jpdiff+b))
%         for(idiff=-a:a)
%             for(jdiff=-a:a)
                idiff2=idiff-ipdiff;
                jdiff2=jdiff-jpdiff;
                if((abs(idiff)+abs(jdiff)<=a) && (abs(idiff2)+abs(jdiff2)<=b))
                    AprodBpmxA(ipdiff+amaxABp+1,jpdiff+amaxABp+1,AprodBpnb(ipdiff+amaxABp+1,jpdiff+amaxABp+1)+(1:3))=invindlist(idiff+amaxABp+1,jdiff+amaxABp+1,1:3);
                    AprodBpmxB(ipdiff+amaxABp+1,jpdiff+amaxABp+1,AprodBpnb(ipdiff+amaxABp+1,jpdiff+amaxABp+1)+(1:3))=invindlist(idiff2+amaxABp+1,jdiff2+amaxABp+1,1:3);
                    AprodBpnb(ipdiff+amaxABp+1,jpdiff+amaxABp+1)=AprodBpnb(ipdiff+amaxABp+1,jpdiff+amaxABp+1)+3;
                end
            end
        end
    
    end
end

end
function ABp=AprodBp(A,B,d,amaxABp,AprodBpmxA,AprodBpmxB,AprodBpnb)
global invindlist;
global indmx;

d2=d^2;d23=d2*3;
% indmaxA=3*(1+amaxA*(amaxA+1)*2);
% indmaxB=3*(1+amaxB*(amaxB+1)*2);
indmaxABp=3*(1+amaxABp*(amaxABp+1)*2);
indmaxB=size(B,2);

ABp=zeros(d23,indmaxABp);
for ipdiff=-amaxABp:amaxABp
    for jpdiff=(-amaxABp+abs(ipdiff)):(amaxABp-abs(ipdiff))
        ipd=ipdiff+amaxABp+1;
        jpd=jpdiff+amaxABp+1;
        nb=AprodBpnb(ipd,jpd);
        for(uvhp=1:3)
            Bp=B(indmx((uvhp-1)*d2+1:uvhp*d2,invindlist(ipdiff+amaxABp+1,jpdiff+amaxABp+1,uvhp)),1:indmaxB);
            for(uvh=1:3)
            %ABp((uvh-1)*d2+1:(uvh)*d2,invindlist(ipdiff+amaxABp+1,jpdiff+amaxABp+1,uvhp))=sum(A((uvh-1)*d2+1:(uvh)*d2,AprodBpmxA(ipd,jpd,1:nb)).*B((uvhp-1)*d2+1:(uvhp)*d2,AprodBpmxB(ipd,jpd,1:nb)),2);
            %ABp((uvh-1)*d2+1:(uvh)*d2,invindlist(ipdiff+amaxABp+1,jpdiff+amaxABp+1,uvhp))=sum(A((uvh-1)*d2+1:(uvh)*d2,AprodBpmxA(ipd,jpd,1:nb)).*Bp((uvhp-1)*d2+1:(uvhp)*d2,AprodBpmxB(ipd,jpd,1:nb)),2);
            ABp((uvh-1)*d2+1:(uvh)*d2,invindlist(ipdiff+amaxABp+1,jpdiff+amaxABp+1,uvhp))=sum(A((uvh-1)*d2+1:(uvh)*d2,AprodBpmxA(ipd,jpd,1:nb)).*Bp(1:d2,AprodBpmxB(ipd,jpd,1:nb)),2);
            end
        end
    end
end


end
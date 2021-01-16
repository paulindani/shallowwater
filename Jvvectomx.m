function [JvmxI,JvmxJ,JvmxV]=Jvvectomx(HdiffT,d,k,N,rx,ry,i,j,uvh)
assert(isa(d, 'double'));
assert(isa(k, 'double'));
assert(isa(N, 'double'));
assert(isa(i, 'double'));
assert(isa(j, 'double'));
assert(isa(rx, 'double'));
assert(isa(ry, 'double'));
assert(isa(uvh, 'double'));
assert(isa(HdiffT, 'double'));
assert(all(size(HdiffT)<=[1001,3000000]));


d2=d^2;d22=2*d2;d23=d2*3;
JvmxI=zeros(k+1,d23);
JvmxJ=zeros(k+1,d23);
JvmxV=zeros(k+1,d23);
% x=modd(1+i,d);
% y=modd(1+j,d);
% xy=(x-1)*d+y;
% 
% grsum(xy)=sum(gsmdiff,2);
N21=2*N+1;

it=0;
for(l=0:rx-1)
    for(m=0:ry-1)
        
        x=modd(1+i+l*(N21),d);
        y=modd(1+j+m*(N21),d);
        xy=(x-1)*d+y+uvh*d2;

        for(xdiff=-N:N)
            for(ydiff=-N:N)
                for(uvhp=0:2)
                it=it+1;
                xp=modd(x+xdiff,d);
                yp=modd(y+ydiff,d);
                xyp=(xp-1)*d+yp+uvhp*d2;

                JvmxI(1:k+1,it)=xyp;
                JvmxJ(1:k+1,it)=xy;
                JvmxV(1:k+1,it)=HdiffT(1:k+1,xyp);
                end
            end
        end
        
    end
end
end
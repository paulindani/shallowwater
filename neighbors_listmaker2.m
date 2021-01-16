function NM=neighbors_listmaker2(d,N)
assert(isa(d, 'double'));
assert(isa(N, 'double'));
d2=d^2;
N21=2*N+1;
NM=zeros(d2,N21^2*3);

for(i=1:d)
    for(j=1:d)
        ij=(i-1)*d+j;
        ijit=0;
        for(idiff=-N:N)
            for(jdiff=-N:N)
                i2=modd(i+idiff,d);
                j2=modd(j+jdiff,d);
                ij2=(i2-1)*d+j2;
                ijit=ijit+1;
                NM(ij,ijit)=ij2;
            end
        end
    end
end
NM(1:d2,N21^2+1:2*N21^2)=NM(1:d2,1:N21^2)+d2;
NM(1:d2,2*N21^2+1:3*N21^2)=NM(1:d2,1:N21^2)+2*d2;
end
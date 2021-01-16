function Dvec=Dvecgen(d,N,rx,ry,i,j)
assert(isa(d, 'double'));
assert(isa(N, 'double'));
assert(isa(i, 'double'));
assert(isa(j, 'double'));
assert(isa(rx, 'double'));
assert(isa(ry, 'double'));

Dvec=zeros(1,d^2);
N21=2*N+1;
for(l=1:rx)
    for(m=1:ry)
        x=modd(1+i+l*(N21),d);
        y=modd(1+j+m*(N21),d);
        Dvec((x-1)*d+y)=1;
    end
end

end
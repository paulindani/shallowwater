function [Ivec,Jvec,Vvec]=neighbors_listmaker0(d)
assert(isa(d, 'double'));
d2=d^2;

Ivec=zeros(3*d2,1);
Jvec=zeros(3*d2,1);
Vvec=ones(3*d2,1);
it=1;
for(i=1:d)
    for(j=1:d)
        ij=(i-1)*d+j;
        d2ij=d2+ij;
        td2ij=2*d2+ij;
        Ivec(it)=ij;Jvec(it)=ij;it=it+1;
        Ivec(it)=ij;Jvec(it)=d2ij;it=it+1;
        Ivec(it)=ij;Jvec(it)=td2ij;it=it+1;
    end
end

end
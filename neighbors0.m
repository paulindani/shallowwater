function M=neighbors0(d)
d2=d^2;d23=3*d2;

[Ivec,Jvec,Vvec]=neighbors_listmaker0_mex(d);
M=sparse(Ivec,Jvec,Vvec,d23,d23,3*d23);

M(d2+1:2*d2,1:d23)=M(1:d2,1:d23);
M(2*d2+1:d23,1:d23)=M(1:d2,1:d23);
end
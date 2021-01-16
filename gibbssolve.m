function X=gibbssolve(M,b,X0,sims,d,nG,nN,MN,MG,MNG)

% assert(isa(X0, 'double'));
%   assert(isa(RM, 'double'));
%   assert(isa(b, 'double'));
%   assert(isa(sims, 'double'));
%   assert(isa(d, 'double'));
%   assert(isa(nG, 'double'));
%   assert(isa(nN, 'double'));
%   assert(isa(MG, 'double'));
%   assert(isa(MNG, 'double'));
% 
%   assert(all(size(X0)<=[2000000,1]));
%   assert(all(size(RM)<=[2000000,1000]));
%   assert(all(size(b)<=[2000000,1]));
%   assert(all(size(MG)<=[2000000,1000]));
%   assert(all(size(MNG)<=[2000000,1000]));



d2=d^2;d23=3*d2;
nN2=(2*nN+1);nN22=nN2^2;
RM=zeros(d23,nN22*3);

for(uvh=1:3)
for(it=1:d2)
    RM(it+(uvh-1)*d2,1:3*nN22)=M(it+(uvh-1)*d2,MN(it,1:3*nN22));
end
end
X=gibbssolve1_mex(X0,RM,b,sims,d,nG,nN,MG,MNG);
end

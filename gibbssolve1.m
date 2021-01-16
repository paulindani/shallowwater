function X=gibbssolve1(X0,RM,b,sims,d,nG,nN,MG,MNG)

assert(isa(X0, 'double'));
  assert(isa(RM, 'double'));
  assert(isa(b, 'double'));
  assert(isa(sims, 'double'));
  assert(isa(d, 'double'));
  assert(isa(nG, 'double'));
  assert(isa(nN, 'double'));
  assert(isa(MG, 'double'));
  assert(isa(MNG, 'double'));

  assert(all(size(X0)<=[2000000,1]));
  assert(all(size(RM)<=[2000000,1000]));
  assert(all(size(b)<=[2000000,1]));
  assert(all(size(MG)<=[2000000,1000]));
  assert(all(size(MNG)<=[2000000,1000]));

d2=d^2;d23=d2*3;
nN2=(2*nN+1);nN22=nN2^2;
nG2=(2*nG+1);nG22=nG2^2; 
nNG2=2*(nN+nG)+1;nNG22=(2*(nN+nG)+1)^2;

X=X0;
for(it=1:sims)
for(i=1:nG:d)
for(j=1:nG:d)
    i=floor(rand*d)+1;
    j=floor(rand*d)+1;
    ij=(i-1)*d+j;
    nGij=MG(ij,1:3*nG22);
    nNGij=MNG(ij,1:3*nNG22);
    MnGnNG=zeros(3*nG22,3*nNG22);
    MnGnG=zeros(3*nG22,3*nG22);
    
    for(idiff=-nG:nG)
        for(jdiff=-nG:nG)
            C=zeros(nNG2);
            for(uvh=1:3)
                ijdiff=(idiff+nG)*nG2+jdiff+nG+1+(uvh-1)*nG22;
                i2=modd(i+idiff,d);
                j2=modd(j+jdiff,d);
                ij2=(i2-1)*d+j2+(uvh-1)*d2;
                C(nG+jdiff+1:nG+jdiff+1+2*nN,nG+idiff+1:nG+idiff+1+2*nN)=...
                    reshape(RM(ij2,1:nN22),[nN2,nN2]);
                R=reshape(C,[1,nNG22]);
                MnGnNG(ijdiff,1:nNG22)=R(1,1:nNG22);%reshape(C,[1,nNG22]);
                R2=reshape(C(nN+1:nN+2*nG+1,nN+1:nN+2*nG+1),[1,nG22]);
                MnGnG(ijdiff,1:nG22)=R2(1,1:nG22);%reshape(C(nN+1:nN+2*nG+1,nN+1:nN+2*nG+1),[1,nG22]);
                C(nG+jdiff+1:nG+jdiff+1+2*nN,nG+idiff+1:nG+idiff+1+2*nN)=...
                    reshape(RM(ij2,nN22+1:2*nN22),[nN2,nN2]);
                R=reshape(C,[1,nNG22]);
                MnGnNG(ijdiff,nNG22+1:2*nNG22)=R(1,1:nNG22);%reshape(C,[1,nNG22]);
                R2=reshape(C(nN+1:nN+2*nG+1,nN+1:nN+2*nG+1),[1,nG22]);
                MnGnG(ijdiff,nG22+1:2*nG22)=R2(1,1:nG22);%reshape(C(nN+1:nN+2*nG+1,nN+1:nN+2*nG+1),[1,nG22]);
                C(nG+jdiff+1:nG+jdiff+1+2*nN,nG+idiff+1:nG+idiff+1+2*nN)=...
                    reshape(RM(ij2,2*nN22+1:3*nN22),[nN2,nN2]);
                R=reshape(C,[1,nNG22]);
                MnGnNG(ijdiff,2*nNG22+1:3*nNG22)=R(1,1:nNG22);%reshape(C,[1,nNG22]);
                R2=reshape(C(nN+1:nN+2*nG+1,nN+1:nN+2*nG+1),[1,nG22]);
                MnGnG(ijdiff,2*nG22+1:3*nG22)=R2(1,1:nG22);%reshape(C(nN+1:nN+2*nG+1,nN+1:nN+2*nG+1),[1,nG22]);
            end
        end
    end
   
    bnnew=b(nGij,1)-MnGnNG*X(nNGij,1);
    X(nGij,1)=X(nGij,1)+MnGnG\bnnew;
end
end
end
end
function [CI,CJ,CV]=localizedcov2(d,Np,neighborsMarr,Vn,ln,neighborstrunc)
  assert(isa(d, 'double'));
  assert(isa(Np, 'double'));
  assert(isa(ln, 'double'));
  assert(isa(neighborsMarr, 'double'));
  assert(isa(Vn, 'double'));
  assert(isa(neighborstrunc, 'double'));
  assert(all(size(neighborsMarr)<=[400000,1000]));
  assert(all(size(Vn)<=[20000,400000]));

  d2=d^2;d23=d2*3;
  CI=zeros(d23*ln,1);
  CJ=zeros(d23*ln,1);
  CV=zeros(d23*ln,1);
  
  it=0;
    for(ind=1:d23)
        for(itind2=1:ln)
            it=it+1;
            ind2=neighborsMarr(ind,itind2);          
            indd2=modd(ind,d2);
            indx=modd(indd2,d);
            indy=(indd2-indx)/d+1;
            
            ind2d2=modd(ind2,d2);
            ind2x=modd(ind2d2,d);
            ind2y=(ind2d2-ind2x)/d+1;

            diffx=abs(indx-ind2x);
            diffx=min(diffx,d-diffx);
            
            diffy=abs(indy-ind2y);
            diffy=min(diffy,d-diffy);
            diff=diffx+diffy;
            
            CI(it)=ind;
            CJ(it)=ind2;
            CV(it)=Vn(1:Np,ind)'*Vn(1:Np,ind2)*rho(2*diff/neighborstrunc);           
        end
    end
end
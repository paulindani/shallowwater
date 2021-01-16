function [CI,CJ,CV]=localizedcov(d,Np,neighborsMarr,Vn,ln)
  assert(isa(d, 'double'));
  assert(isa(Np, 'double'));
  assert(isa(ln, 'double'));
  assert(isa(neighborsMarr, 'double'));
  assert(isa(Vn, 'double'));
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
            CI(it)=ind;
            CJ(it)=ind2;
            CV(it)=Vn(1:Np,ind)'*Vn(1:Np,ind2);           
        end
    end
end
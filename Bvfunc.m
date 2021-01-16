function Bv=Bvfunc(ur,d)
d23=3*d^2;


Bvsp=Bvfuncsp(ur,d);
Bv=sparse(Bvsp(:,2),Bvsp(:,1),Bvsp(:,3),d23,d23);

% Bv=sparse(d23,d23);
% 
% for(it=1:32*d2)
% Bv(Blist(it,2),Blist(it,1))=Bv(Blist(it,2),Blist(it,1))+Blist(it,4)*ur(Blist(it,3));
% end

end
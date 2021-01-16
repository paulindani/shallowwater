function spBv=Bvfuncsp(ur,d)
global Blist4;
global Blist2;
global Blist1;
%global d;
assert(isa(d,'double'));
assert(isa(ur,'double'));
assert(all(size(ur)<=[1,4000000]));

d2=d^2;

spBv=zeros(23*d2,3);

spBv(1:d2,1:2)=Blist4(1:4:(4*d2-3),1:2);
spBv(1:d2,3)=Blist4(1:4:(4*d2-3),4).*ur(Blist4(1:4:(4*d2-3),3))'+Blist4(2:4:(4*d2-2),4).*ur(Blist4(2:4:(4*d2-2),3))'+Blist4(3:4:(4*d2-1),4).*ur(Blist4(3:4:(4*d2-1),3))'+Blist4(4:4:(4*d2),4).*ur(Blist4(4:4:(4*d2),3))';


spBv(d2+1:7*d2,1:2)=Blist2(1:2:(12*d2-1),1:2);
spBv(d2+1:7*d2,3)=Blist2(1:2:(12*d2-1),4).*ur(Blist2(1:2:(12*d2-1),3))'+Blist2(2:2:(12*d2),4).*ur(Blist2(2:2:(12*d2),3))';

spBv(7*d2+1:23*d2,1:2)=Blist1(1:16*d2,1:2);
spBv(7*d2+1:23*d2,3)=Blist1(1:16*d2,4).*ur(Blist1(1:16*d2,3))';

% for(it=1:2*d2)
%     itm14p1=(it-1)*4+1;
%     spBv(it,1:2)=Blist4(itm14p1,1:2);
%     spBv(it,3)=Blist4(itm14p1,4)*ur(Blist4(itm14p1,3))+Blist4(itm14p1+1,4)*ur(Blist4(itm14p1+1,3))+Blist4(itm14p1+2,4)*ur(Blist4(itm14p1+2,3))+Blist4(itm14p1+3,4)*ur(Blist4(itm14p1+3,3));
% end
% 
% for(it=1:4*d2)
%     itm12p1=(it-1)*2+1;
%     spBv(2*d2+it,1:2)=Blist2(itm12p1,1:2);
%     spBv(2*d2+it,3)=Blist2(itm12p1,4)*ur(Blist2(itm12p1,3))+Blist2(itm12p1+1,4)*ur(Blist2(itm12p1+1,3));
% end

% for(it=1:16*d2)
%     spBv(6*d2+it,1:2)=Blist1(it,1:2);
%     spBv(6*d2+it,3)=Blist1(it,4)*ur(Blist1(it,3));
% end



% for(it=1:2*d2)
%     itm14p1=(it-1)*4+1;
%     spBv(it,1:2)=Blist4(itm14p1,1:2);
%     spBv(it,3)=Blist4(itm14p1,4)*ur(Blist4(itm14p1,3))+Blist4(itm14p1+1,4)*ur(Blist4(itm14p1+1,3))+Blist4(itm14p1+2,4)*ur(Blist4(itm14p1+2,3))+Blist4(itm14p1+3,4)*ur(Blist4(itm14p1+3,3));
% end
% 
% for(it=1:4*d2)
%     itm12p1=(it-1)*2+1;
%     spBv(2*d2+it,1:2)=Blist2(itm12p1,1:2);
%     spBv(2*d2+it,3)=Blist2(itm12p1,4)*ur(Blist2(itm12p1,3))+Blist2(itm12p1+1,4)*ur(Blist2(itm12p1+1,3));
% end
% 
% for(it=1:16*d2)
%     spBv(6*d2+it,1:2)=Blist1(it,1:2);
%     spBv(6*d2+it,3)=Blist1(it,4)*ur(Blist1(it,3));
% end


end
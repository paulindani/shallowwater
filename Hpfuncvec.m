 function v=Hpfuncvec(vobs)
%   v=zeros(1,d);
%   for i = 1:div(d,6)
%     v(6*(i-1)+1:6*(i-1)+3)=vobs(3*(i-1)+1:3*(i-1)+3);
%   end
    %v=vobs;
    global d;
    global obsuvprec;

    v=obsuvprec.*vobs;
%    d23=3*d^2;
%     v=zeros(1,d23);
%     v(2*d^2+1:d23)=vobs;
 end
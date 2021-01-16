function [gsmval,grad,grad2]=gradhessgsm(imax,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)
%global dobs;
global FactorialArray;
global BinomialArray;

d23=3*d^2;
dobs=d23;

gsmval=0;
grad=zeros(1,d23);
grad2=sparse(d23,d23);
%grad2=spalloc(d23,d23,4*(imax+1)^2*d23);
Deltaobs=zeros(k+1,dobs);

vT = zeros(k+1,d23);
vT(1, 1:d23) = v(1:d23);

vder=dercalc(imax, v,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray);
vderJv=derJvcalc(imax, v, vder,d,H,Hip1,Him1,Hjp1,Hjm1);

%vderJv2=derJv2calc(imax, v, d, F, BinomialArray,vder,vderJv);


  TFarr=zeros(k,imax+1);
  
  for it = 1:k
      for j = 1:imax+1
        TFarr(it,j)=((T(it + 1))^(j-1))/FactorialArray(j);
      end
  end
  
  vT(2:k+1,1:d23)=TFarr*vder;

% for it = 1:k
%   vT(it+1,1:d23)=v(1:d23);
%   for j = 1:imax
%      vT(it+1,1:d23)=vT(it+1,1:d23)+vder(j + 1,1:d23)*(T(it + 1))^j/FactorialArray(j + 1);
%   end
% end

PhiT=Hfuncmx(vT,k,d);

Deltaobs=Y-PhiT;

gsmval=sum(sum(Deltaobs.^2))/2;

Deltapow=zeros(imax+1,dobs);
HpDeltapow=zeros(imax+1,d23);

imax21=2*imax+1;
Tpow=zeros(imax21,k+1);
Tpow(1,1:k+1)=1;
for i = 1:2*imax
Tpow(i+1,1:k+1)=T.^i;%.*Tpow(i,1:k+1);
end

for i = 0:imax
  for j = 1:k+1
    Deltapow(i+1,1:dobs)=Deltapow(i+1,1:dobs)+Tpow(i+1,j)*Deltaobs(j,1:dobs);
  end
end


for i = 0:imax
  HpDeltapow(i+1,1:d23)=Hpfuncvec(Deltapow(i+1,1:dobs));
    %grad=grad-circshift(vderJv(i+1,l,1:d23).*HpDeltapow(i+1,1:d23),-(2*imax+1-l))/FactorialArray(i+1);
    grad=grad-HpDeltapow(i+1,1:d23)*vderJv{i+1}/FactorialArray(i+1);
end


%return;




% for i = 0:imax
%   for i6 = 0:div(d,6)-1
%     for iplus = 1:3
%       i1=i6*6+iplus;i1obs=3*i6+iplus;
%     for i2mi1 = -2*imax:2*imax
%       for i3mi1 = -2*imax:2*imax
%         grad2(i3mi1-i2mi1+imax41,modd(i1+i2mi1,d))-=vderJv2(i+1,i2mi1+imax21,i3mi1+imax21,i1)*Deltapow(i+1,i1obs)/FactorialArray(i+1);
%       end
%     end
%   end
% end
% end


Tpowsum=zeros(imax21,1);
for i = 0:2*imax
  Tpowsum(i+1)=sum(Tpow(i+1,1:k+1));
end



% for i = 0:imax
%   for j = 0:imax
%     %grad2=grad2+Tpowsum(i+j+1)*HApBfunc(vderJv{i+1}, vderJv{j+1})/(FactorialArray(i+1)*FactorialArray(j+1));
%     grad2=grad2+Tpowsum(i+j+1)*(vderJv{i+1})'*vderJv{j+1}/(FactorialArray(i+1)*FactorialArray(j+1));
%     %grad2=grad2+(Tpowsum(i+j+1)/(FactorialArray(i+1)*FactorialArray(j+1)))*spmxprod((vderJv{i+1})',vderJv{j+1},2*imax);
%   end
% end

%cell(imax+1,1);
%for(it=1:imax+1)
%sumarr{it}=spalloc(d23,d23,4*(imax+1)^2*d23);
%end

for i = 0:imax
  %sumarr = sparse(d23,d23);%,4*(imax+1)^2*d23);
  sumarr=spalloc(d23,d23,4*(imax+1)^2*d23);
  for j = 0:imax
    %grad2=grad2+Tpowsum(i+j+1)*HApBfunc(vderJv{i+1}, vderJv{j+1})/(FactorialArray(i+1)*FactorialArray(j+1));
    sumarr=sumarr+(Tpowsum(i+j+1)/(FactorialArray(i+1)*FactorialArray(j+1)))*vderJv{j+1};
    %grad2=grad2+(Tpowsum(i+j+1)/(FactorialArray(i+1)*FactorialArray(j+1)))*spmxprod((vderJv{i+1})',vderJv{j+1},2*imax);
  end
  %grad2=grad2+(vderJv{i+1})'*sumarr;
  grad2=grad2+HApBfunc(vderJv{i+1},sumarr);
end

%sDeltaobs=sum(Deltaobs,1);
% sDeltaobsd23=zeros(d23,1);
% Tdobs=T'*ones(1,dobs);
% sDeltaobsd23(1:d23)=sum(Deltaobs.*Tdobs,1);


%grad2=grad2+2*Bvmiddlefunc(sDeltaobsd23');



% for i = 0:imax-1
%   for j = i+1:imax
%     grad2=grad2+Tpowsum(i+j+1)*(vderJv{i+1})'*vderJv{j+1}/(FactorialArray(i+1)*FactorialArray(j+1));
%   end
% end
% grad2=grad2+grad2';
% for i = 0:imax
%     j=i;
%     grad2=grad2+Tpowsum(i+j+1)*(vderJv{i+1})'*vderJv{j+1}/(FactorialArray(i+1)*FactorialArray(j+1));
% end




%return gsmval,grad,grad2;

end



 
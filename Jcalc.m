function J=Jcalc(imax,v,Y,T,k,d,H,Hip1,Him1,Hjp1,Hjm1)
%global dobs;
global FactorialArray;
global BinomialArray;

d23=3*d^2;
dobs=d23;

gsmval=0;
grad=zeros(1,d23);
%grad2=sparse(d23,d23);
J=spalloc(d23,d23,15*(imax)^2*d23);
Deltaobs=zeros(k+1,dobs);


vder=dercalc(imax, v,d,H,Hip1,Him1,Hjp1,Hjm1);
vderJv=derJvcalc(imax, vder,d);




  TFarr=zeros(k,imax+1);
  
  for it = 1:k
      for j = 1:imax+1
        TFarr(it,j)=((T(it + 1))^(j-1))/FactorialArray(j);
      end
  end
  
  for(it=1:imax+1)
      J=J+TFarr(k,it)*vderJv{it};
  end
  
 

end



 
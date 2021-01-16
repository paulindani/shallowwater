function uT=onestep(u0, T, k,dermax,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray,FactorialArray)
  %global BinomialArray;
  %global FactorialArray;
  %global d;
  
  
  assert(isa(dermax, 'double'));
  assert(isa(u0, 'double'));
  assert(isa(k, 'double'));
  assert(isa(d, 'double'));
  assert(isa(T, 'double'));
  assert(isa(H, 'double'));
  assert(isa(Hip1, 'double'));
  assert(isa(Him1, 'double'));
  assert(isa(Hjp1, 'double'));
  assert(isa(Hjm1, 'double'));
  assert(isa(BinomialArray, 'double'));
  assert(isa(FactorialArray, 'double'));

  assert(all(size(u0)<=[1,4000000]));
  assert(all(size(H)<=[1,4000000]));
  assert(all(size(Hip1)<=[1,4000000]));
  assert(all(size(Him1)<=[1,4000000]));
  assert(all(size(Hjp1)<=[1,4000000]));
  assert(all(size(Hjm1)<=[1,4000000]));
  assert(all(size(BinomialArray)<=[51,51]));
  assert(all(size(FactorialArray)<=[1,51]));

  assert(all(size(T)<=[1,10001]));
  
  d23=3*d^2;
  uder =dercalc(dermax, u0,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray);

  uT = zeros(k+1,d23);
  uT(1, 1:d23) = u0(1:d23);

  TFarr=zeros(k,dermax+1);
  
  for it = 1:k
      for j = 1:dermax+1
        TFarr(it,j)=((T(it + 1))^(j-1))/FactorialArray(j);
      end
  end
  uT(2:k+1,1:d23)=TFarr*uder;
 
end

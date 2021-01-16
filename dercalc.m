function uder=dercalc(imax, u0,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray)
  %global BinomialArray;
  assert(isa(imax, 'double'));
  assert(isa(d, 'double'));
  assert(isa(u0, 'double'));
  assert(all(size(u0)<=[1,4000000]));
  assert(isa(BinomialArray, 'double'));
  assert(all(size(BinomialArray)<=[51,51]));

  assert(isa(H, 'double'));
  assert(isa(Hip1, 'double'));
  assert(isa(Him1, 'double'));
  assert(isa(Hjp1, 'double'));
  assert(isa(Hjm1, 'double'));

  assert(all(size(H)<=[1,4000000]));
  assert(all(size(Hip1)<=[1,4000000]));
  assert(all(size(Him1)<=[1,4000000]));
  assert(all(size(Hjp1)<=[1,4000000]));
  assert(all(size(Hjm1)<=[1,4000000]));

 
  d23=3*d^2;
  uder = zeros(imax+1,d23);
  uder(1,1:d23)=u0(1,1:d23);
  uder(2,1:d23) = -Afunc(u0,d,H,Hip1,Him1,Hjp1,Hjm1) - Bfunc(u0,u0,d);
%  uder(3,1:d23) = -Afunc(uder(2,1:d23)) - Bfunc(uder(2,1:d23),u0)-Bfunc(u0,uder(2,1:d23));
  for it=2:imax
    uderit=-Afunc(uder(it,1:d23),d,H,Hip1,Him1,Hjp1,Hjm1);
%     for j =0:it-1
%         uderit=uderit-BinomialArray(it,j+1)*Bfunc(uder(it - j,1:d23),uder(j + 1,1:d23));
%     end
    if(mod(it-1,2)==1)
      for j =0:div(it-2,2)
        uderit=uderit-2*BinomialArray(it,j + 1)*Bfunc(uder(it - j,1:d23),uder(j + 1,1:d23),d);
      end
    else
      for j= 0:div(it-3,2)
        uderit=uderit-2*BinomialArray(it,j + 1)*Bfunc(uder(it - j,1:d23),uder(j + 1,1:d23),d);
      end
      j=div(it-1,2);
      uderit=uderit-BinomialArray(it,j + 1)*Bfunc(uder(it - j,1:d23),uder(j + 1,1:d23),d);
    end
    uder(it + 1, 1:d23) = uderit;
  end
end


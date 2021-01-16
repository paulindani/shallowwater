function BinomialArray=BinomialArrayFunc(m)
  global BinomialArray;
  BinomialArray=zeros(m+1,m+1);
  for i1= 0:m
    for i2= 0:i1
      BinomialArray(i1+1,i2+1)=nchoosek(i1,i2);
    end
  end
end
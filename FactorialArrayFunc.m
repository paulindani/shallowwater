function FactorialArrayFunc(m)
   global FactorialArray;
   FactorialArray=zeros(1,m+1);
   for i= 0:m
       FactorialArray(i+1)=factorial(i);
   end
end


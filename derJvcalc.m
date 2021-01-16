 function uderJv=derJvcalc(imax, uder,d)
  global BinomialArray;
  %global d;
  %global neighborsM;
  global Amx;
  
  d23=d^2*3;
  uderJv = cell(imax+1,1);
  Bdv = cell(imax+1,1);
  for(i=1:imax+1)
      uderJv{i}=spalloc(d23,d23,20*d23);
      Bdv{i}=Bvfunc(uder(i,1:d23),d);
  end
  uderJv{1}=speye(d23);
  
  %zeros(imax+1,imax41,d);
  %uderJv{2}=-Afuncmx(uderJv{1})-Bvfunc(u0);
  for it = 2:imax+1
     uderJvit=-Amx*uderJv{it-1};%Afuncmx(uderJv{it-1});
     for j = 0:it-2
        %uderJvit=uderJvit-(2*BinomialArray(it-1,j+1))*BfuncJvdv(uderJv{j+1},uder(it-1-j,1:d23));
        uderJvit=uderJvit-(2*BinomialArray(it-1,j+1))*Bdv{it-1-j}*uderJv{j+1};
     end
     uderJv{it}=uderJvit;
     %uderJv{it} = neighborsM.*uderJvit;
  end
  end
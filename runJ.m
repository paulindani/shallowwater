function [uT,Jarr,Jinvarr]=runJ(imaxu,imaxJ, u0, T, k,d,H,Hip1,Him1,Hjp1,Hjm1)
global FactorialArray;
global BinomialArray;
global interpol;
global neighborsM;


d23=3*d^2;
uT=zeros(k+1,d23);
uT(1,1:d23)=u0;
Jarr=cell(k,1);
if(nargout==3)Jinvarr=cell(k,1);end

  TFarr=zeros(1,max(imaxu,imaxJ)+1);
  
  h=T(2)-T(1);
  for j = 1:max(imaxu,imaxJ)+1
    TFarr(j)=(h)^(j-1)/FactorialArray(j);
  end
  
  
  if(nargout==2)
      for (l=1:k)
          uder =dercalc_mex(imaxu, uT(l,1:d23),d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray);
          uT(l+1,1:d23)=TFarr(1,1:imaxu+1)*uder;
          if(l==k|| l==1|| l==2 || mod(l,interpol)==0)
          uderJ=derJvcalc(imaxJ, uder,d);  
          Jarr{l}=spalloc(d23,d23,39*d23);
          
          Jarrl=spalloc(d23,d23,3*(1+2*imaxJ*(imaxJ+1))*d23);
          for (it=1:imaxJ+1)
              %Jarr{l}=Jarr{l}+TFarr(it)*uderJ{it};
              Jarrl=Jarrl+TFarr(it)*uderJ{it};
          end
          Jarr{l}=neighborsM.*Jarrl;
          %Jarr{l}=Jarrl;%neighborsM.*Jarrl;
          end
      end
  
  else   
      for(l=1:k)
          uder =dercalc_mex(imaxu, uT(l,1:d23),d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray);
          uT(l+1,1:d23)=TFarr(1,1:imaxu+1)*uder;
          uderJ=derJvcalc(imaxJ, uder,d);
          if(l==k|| l==1|| l==2|| mod(l,interpol)==0)
          Jarr{l}=spalloc(d23,d23,39*d23);
          if(l>1)Jinvarr{l-1}=spalloc(d23,d23,39*d23);end
          Jarrl=spalloc(d23,d23,3*(1+2*imaxJ*(imaxJ+1))*d23);
          Jinvarrlm1=spalloc(d23,d23,3*(1+2*imaxJ*(imaxJ+1))*d23);
          for(it=1:imaxJ+1)
%               Jarr{l}=Jarr{l}+TFarr(it)*uderJ{it};
%               if(l>1)Jinvarr{l-1}=Jinvarr{l-1}+TFarr(it)*(-1)^(it+1)*uderJ{it};end
              Jarrl=Jarrl+TFarr(it)*uderJ{it};
              if(l>1)Jinvarrlm1=Jinvarrlm1+TFarr(it)*(-1)^(it+1)*uderJ{it};end
          end
          Jarr{l}=neighborsM.*Jarrl;%Jarrl;
          if(l>1)Jinvarr{l-1}=neighborsM.*Jinvarrlm1;%Jinvarrlm1;
          end
%           Jarr{l}=neighborsM.*Jarr{l};
%           if(l>1)Jinvarr{l-1}=neighborsM.*Jinvarr{l-1};end
          %end
      end
      uder =dercalc_mex(imaxu, uT(k+1,1:d23),d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray);
      uderJ=derJvcalc(imaxJ, uder,d);  
      Jinvarrk=spalloc(d23,d23,39*d23);
      for(it=1:imaxJ+1)
          Jinvarrk=Jinvarrk+TFarr(it)*(-1)^(it+1)*uderJ{it};
      end
      Jinvarr{k}=neighborsM.*Jinvarrk;
  end
   


end
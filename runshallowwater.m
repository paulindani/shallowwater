function uT=runshallowwater(dermax, nbblocks, u0, T, k,d,H,Hip1,Him1,Hjp1,Hjm1)%,BArray,FArray,ddd)
global BinomialArray;
global FactorialArray;
%global d;
%global H;
%d2=d^2;
% BinomialArray=BArray;
% FactorialArray=FArray;
% d=ddd;

% assert(isa(dermax, 'int64'));
% assert(isa(nbblocks, 'int64'));
% assert(isa(k, 'int64'));
% assert(isa(T, 'double'));
% assert(all(size(T)<=[101]));

% 

if(k==0)
uT=u0;
return;
end

if(nbblocks==1)
    uT=onestep_mex(u0, T, k,dermax,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray,FactorialArray);
    return;
end
    
d23=3*d^2;
uT=zeros(k+1,d23);


if(nbblocks<=k)
if(mod(k,nbblocks)==0)
kblockarr=(k/nbblocks)*ones(1,nbblocks);
else
kblockarr=zeros(1,nbblocks);
fk=floor(k/nbblocks);
b=mod(k,fk);
a=nbblocks-b;
kblockarr(1:a)=fk;
kblockarr(a+1:nbblocks)=fk+1;
end
else
kblockarr=ones(1,k);
end

ckblockarr=cumsum(kblockarr);
ckblockarr=[0, ckblockarr];
u0it=u0;
for it =1:length(kblockarr);
uT(ckblockarr(it)+1:ckblockarr(it+1)+1,1:d23)=onestep_mex(u0it, T(1:kblockarr(it)+1), kblockarr(it),dermax,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray,FactorialArray);
u0it=uT(ckblockarr(it+1)+1,1:d23);
end


end


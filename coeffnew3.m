function dercoeffarr=coeffnew3(lmax,k,Tarr,nbblocks,Y,sigmaz,d,H,Hip1,Him1,Hjp1,Hjm1)
% global BinomialArray;
% global FactorialArray;
% global dobs;
global Delta;
% global f;
% global g;
% global h;
% global obsuvprec;
% global spobsuvprec;
% global neighborsM;
%global sigmaz;
global priorprecmx;
global priormean;

imax=10;N=3;delta=1e-2;

d2=d^2; d23=3*d2;
dercoeffarr=zeros(lmax+1,k+1);

optionslsqlin=optimoptions('lsqlin','Algorithm','interior-point','MaxIter',10000,'TolCon',2*eps,'TolFun',2*eps);
%optionslsqlin=optimoptions('lsqlin','MaxIter',1000);
    
 
%     vest=Y(1,1:d23);
%     [gsmval,grad,grad2]=gradhessgsmnew(imax,nbblocks,N,delta,vest,Y(1,1:d23),Tarr(1),0,d,H,Hip1,Him1,Hjp1,Hjm1);
%     vest=vest-(cholmod2(grad2,grad'))'; 
    vest=priormean;
    cholpriorprecmx=chol(priorprecmx)/sigmaz;
    

    numsims=10*k;
    u0arr=zeros(numsims,d23);
    u0arr(1:numsims,1:d23)=ones(numsims,1)*vest+(cholpriorprecmx\randn(d23,numsims))';

    
    M=zeros(numsims,k);
    uTarrsims=zeros(numsims,k+1,d23);
    Ysims=zeros(numsims,k+1,d23);
    udersims=zeros(numsims,lmax+1,d23);
for(it=1:numsims)
    uTarrsims(it,1:k+1,1:d23)=runshallowwater(imax,nbblocks, u0arr(it,1:d23), Tarr, k,d,H,Hip1,Him1,Hjp1,Hjm1);
    Ysims(it,1:k+1,1:d23)=uTarrsims(it,1:k+1,1:d23)+sigmaz*reshape(randn(k+1,d23),[1,k+1,d23]);%.*(ones(k+1,1)*obsuvprec)
    udersims(it,1:lmax+1,1:d23)=dercalc(lmax,u0arr(it,1:d23),d,H,Hip1,Him1,Hjp1,Hjm1);
end

%M=zeros(numsims*d23,(k+1));
% for(ij=1:d23)
%     M((ij-1)*numsims+1:ij*numsims,1:k+1)=reshape(Ysims(1:numsims,1:k+1,ij),[numsims,k+1]);
% end

M=zeros(numsims*d2,(k+1));
for(ij=1:d2)
    M((ij-1)*numsims+1:ij*numsims,1:k+1)=reshape(Ysims(1:numsims,1:k+1,2*d2+ij),[numsims,k+1]);
end


    MpM=M'*M;
%     rMpM=rank(MpM)
    
%     size(MpM)
%     cond(MpM)
    V=zeros(numsims*d2,lmax+1);
    for(derit=1:lmax+1)
    v=reshape(udersims(1:numsims,derit,2*d2+1:d23),[numsims*d2,1])*((Delta)^(derit-1));
%    dercoeffarr(derit,1:(k+1))=((MpM)\(M'*v))/((Delta)^(derit-1));
    %dercoeffarr(derit,1:k+1)=lsqlin(M,v,[],[])/((10*Delta)^(derit-1));
    %dercoeffarr(derit,1:2*(k+1))=lsqlin(M,v,[],[])/((10*Delta)^(derit-1));
    dercoeffarr(derit,1:k+1)=lsqlin(M,v,[],[],[],[],[],[],[],optionslsqlin)/((Delta)^(derit-1));
    end
    
    
    
end    


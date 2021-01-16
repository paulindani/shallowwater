function dercoeffarr=coeffnew4(lmax,k,Tarr,nbblocks,Y,sigmaz,d,H,Hip1,Him1,Hjp1,Hjm1)
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

imax=3;N=3;delta=1e-2;

d2=d^2; d23=3*d2;
dercoeffarr=zeros(lmax+1,k+1);

optionslsqlin=optimoptions('lsqlin','Algorithm','interior-point','MaxIter',10000,'TolCon',2*eps,'TolFun',2*eps);
%optionslsqlin=optimoptions('lsqlin','MaxIter',1000);
    
 
%     vest=Y(1,1:d23);
%     [gsmval,grad,grad2]=gradhessgsmnew(imax,nbblocks,N,delta,vest,Y(1,1:d23),Tarr(1),0,d,H,Hip1,Him1,Hjp1,Hjm1);
%     vest=vest-(cholmod2(grad2,grad'))'; 
    vest=priormean;
    %cholpriorprecmx=chol(priorprecmx)/sigmaz;
    cholprecmx=speye(d23)/sigmaz/10;

    numsims=10*k;
    u0arr=zeros(numsims,d23);
    u0arr(1:numsims,1:d23)=ones(numsims,1)*vest+(cholprecmx\randn(d23,numsims))';

    
    M=zeros(numsims,k);
    uTarrsims=zeros(numsims,k+1);
    Ysims=zeros(numsims,k+1);
    udersims=zeros(numsims,lmax+1);
for(it=1:numsims)
    uTtemp=runshallowwater(imax,nbblocks, u0arr(it,1:d23), Tarr, k,d,H,Hip1,Him1,Hjp1,Hjm1);
    uTarrsims(it,1:k+1)=uTtemp(1:k+1,1);
    Ysims(it,1:k+1)=uTtemp(1:k+1,1)+sigmaz*randn(k+1,1);
    %uTarrsims(it,1:k+1,1:d23)+sigmaz*reshape(randn(k+1,d23),[1,k+1,d23]);%.*(ones(k+1,1)*obsuvprec)
    udersimstemp=dercalc(lmax,u0arr(it,1:d23),d,H,Hip1,Him1,Hjp1,Hjm1);
    udersims(it,1:lmax+1)=udersimstemp(1:lmax+1,1);
end

%M=zeros(numsims*d23,(k+1));
% for(ij=1:d23)
%     M((ij-1)*numsims+1:ij*numsims,1:k+1)=reshape(Ysims(1:numsims,1:k+1,ij),[numsims,k+1]);
% end

M=Ysims;


    MpM=M'*M;
%     rMpM=rank(MpM)
    
%     size(MpM)
%     cond(MpM)
    V=zeros(numsims*d2,lmax+1);
    for(derit=1:lmax+1)
    v=udersims(1:numsims,derit)*((Delta)^(derit-1));
%    dercoeffarr(derit,1:(k+1))=((MpM)\(M'*v))/((Delta)^(derit-1));
    %dercoeffarr(derit,1:k+1)=lsqlin(M,v,[],[])/((10*Delta)^(derit-1));
    %dercoeffarr(derit,1:2*(k+1))=lsqlin(M,v,[],[])/((10*Delta)^(derit-1));
    dercoeffarr(derit,1:k+1)=lsqlin(M,v,[],[],[],[],[],[],[],optionslsqlin)/((Delta)^(derit-1));
    end
    
    
    
end    


function ret=shallowwater_hybrid4dvar_climat(filename,params)
kass=params.kass;
nbblocksass=params.nbblocksass;
interpolin=params.interpolin;
neighborstrunc=params.neighborstrunc;
loc=params.loc;
pcgsteps=params.pcgsteps;
Nparticles=params.Nparticles;
inflation=params.inflation;
hybridc=params.hybridc;

tic();
global sigmaz;
global d;
load(filename);
Tassarr=Tarr(1:kass+1);
Tass=Tassarr(kass+1);


global kcall;
global dcall;
global Jarrcall;

global lastJ;
global lastJinv;
global itJ;
global nbprevJstored;
global mJcall;

global prioruv;
global priorh;
global priorprecmx;
global C;
global priormean;

global Blist4;
global Blist2;
global Blist1;
global Amx;
global interpol;
global neighborsM;


prioruv=params.prioruv;
priorh=params.priorh;

%priorprecen4dvar(d);
%priorprecmx0=priorprecmx;
%C0=C;
[Cmx,Cmean]=climatological_cov(filename);
C0=Cmx;
priorprecmx0=inv(Cmx);
priormean=Cmean;

%%%%%%%Initial definitions 
neighborsM=neighbors(neighborstrunc,d);
neighborsMloc=neighbors(loc,d);

neighborlist=find(neighborsMloc(:,1));
ln=length(neighborlist);
neighborsMarr=zeros(d23,ln);
for(ind=1:d23)
neighborsMarr(ind,1:ln)=find(neighborsMloc(:,ind));
end
Gammamx=sigmaz^2*speye(d23);

imax1=4;
imaxsim = 5;
imaxJ=4;
mJ=0;
mJcall=0;
interpol=interpolin;

    



Amx=genAmx(d,H);
[Blist1,Blist2,Blist4]=genBlist4(d);
lastJ=cell(mJ,1);
lastJinv=cell(mJ,1);
itJ=0;

    






    



    diffuv=zeros(k+1,1);
    diffh=zeros(k+1,1);
    reldiffuv=zeros(k+1,1);
    reldiffh=zeros(k+1,1);  
    nbtot=floor(k/kass);
    err0arr=zeros(nbtot,1);
    
    for(ittot=1:nbtot)
    ittot

    nbprevJstored=0;
    
    if(ittot==1)
        [priormean,priorprecmx]=priorupdateinitialobs(priorprecmx0,priormean,Y(1,1:d23));
        C=inv(priorprecmx);              
        v=priormean;
        %cp=condest(priorprecmx)
        varr=zeros(Nparticles,d23);
        diff=zeros(k+1,1); 
        cholgrad2=chol(C)*sigmaz;
        varr(1:Nparticles,1:d23)=ones(Nparticles,1)*v+(cholgrad2*randn(d23,Nparticles))';    
        mnew=(sum(varr(1:Nparticles,1:d23),1)/Nparticles)';             
        vnew=zeros(Nparticles,d23);
    
        err0arr(1)=sqrt(mean((v(1:2*d2)-u0(1:2*d2)).^2));
        vTarr=runshallowwater(imaxsim,nbblocksass, v, Tassarr,  kass,d,H,Hip1,Him1,Hjp1,Hjm1);
        for(l=1:kass+1)
            diffuv(l)=sqrt(mean((abs(vTarr(l,1:2*d2)-uTarr(l,1:2*d2))).^2));
            diffh(l)=sqrt(mean((abs(vTarr(l,2*d2+1:d23)-uTarr(l,2*d2+1:d23))).^2));
            reldiffuv(l)=diffuv(l)/sqrt(mean((uTarr(l,1:2*d2)).^2));
            reldiffh(l)=diffh(l)/sqrt(mean((uTarr(l,2*d2+1:d23)).^2));
        end
        
    else
        C=C0*hybridc+C*(1-hybridc);
        vfi=pushforward(imaxsim, nbblocksass, v, Tass,d,H,Hip1,Him1,Hjp1,Hjm1);                     
        err0arr(ittot)=sqrt(mean((v(1:2*d2)-uTarr((ittot-1)*kass+1,1:2*d2)).^2));
        vTarr=runshallowwater(imaxsim,nbblocksass, vfi, Tassarr, kass,d,H,Hip1,Him1,Hjp1,Hjm1);
        for(l=1:kass+1)
            diffuv((ittot-1)*kass+l)=sqrt(mean((abs(vTarr(l,1:2*d2)-uTarr((ittot-1)*kass+l,1:2*d2))).^2));
            diffh((ittot-1)*kass+l)=sqrt(mean((abs(vTarr(l,2*d2+1:d23)-uTarr((ittot-1)*kass+l,2*d2+1:d23))).^2));
            reldiffuv((ittot-1)*kass+l)=diffuv((ittot-1)*kass+l)/sqrt(mean((uTarr((ittot-1)*kass+l,1:2*d2)).^2));
            reldiffh((ittot-1)*kass+l)=diffh((ittot-1)*kass+l)/sqrt(mean((uTarr((ittot-1)*kass+l,2*d2+1:d23)).^2));
        end
        v=vfi;
        priormean=vfi;
        
        %recenter particles to vfi;
        vfimestdiff=vfi-mest;
        for(p=1:Nparticles)
        varr(p,1:d23)=varr(p,1:d23)+vfimestdiff;
        end
        
    end

    if(ittot==1)
    gsmvaltrue=gsmJen4dvar(u0,uTarr((ittot-1)*kass+1:ittot*kass+1,1:d23),Y((ittot-1)*kass+1:ittot*kass+1,1:d23))%gsmwithprior(imaxsim,nbblocksass,uTarr((ittot-1)*kass+1,1:d23),Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Tassarr,kass,d,H,Hip1,Him1,Hjp1,Hjm1)
    end
    

    
    
    
    
    stepsize=1;
    Newtonit=0;
    err0=sqrt(mean((v-uTarr((ittot-1)*kass+1,:)).^2))


    
    while stepsize>1e-3 && (Newtonit<5)&& (Newtonit<1|| ittot==1)
    %tic
    
    Newtonit=Newtonit+1
    
    [vTass,Jarr]=runJ(imaxsim,imaxJ, v, Tassarr, kass,d,H,Hip1,Him1,Hjp1,Hjm1);
        
        
        if(ittot==1)
        if(Newtonit==1)gsmval=gsmwithprior(imaxsim,nbblocksass,v,Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Tassarr,kass,d,H,Hip1,Him1,Hjp1,Hjm1);end
        gradnJ=gradJ(v,vTass,Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Jarr,kass,d);   
        else
        if(Newtonit==1)gsmval=gsmwithprioren4dvar(imaxsim,nbblocksass,vn,Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Tassarr,kass,d,H,Hip1,Him1,Hjp1,Hjm1);end
        gradnJ=gradJen4dvar(v,vTass,Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Jarr,kass,d);
        end
            
    %make a pcg step according to the HessJcall that encodes the Hessian
    %based on the information stored in the Jacobian array Jarrcall
    kcall=kass;dcall=d;Jarrcall=Jarr;
    
    if(ittot==1)
        vstep=(pcg(@HessJcall,gradnJ,1e-6,pcgsteps))';
    else
        vstep=(gmres(@HessJcallen4dvar,C*gradnJ,d23,1e-6,pcgsteps))';
    end
 

    vnisgood=0;
    vstep=2*vstep;
    vnit=0;

    while(vnisgood==0 && vnit<=10)
        vnit=vnit+1;
        vstep=vstep/2;
        vn=v-vstep;
        if(ittot==1)
        gsmvn=gsmwithprior(imaxsim,nbblocksass,vn,Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Tassarr,kass,d,H,Hip1,Him1,Hjp1,Hjm1);
        else
        gsmvn=gsmwithprioren4dvar(imaxsim,nbblocksass,vn,Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Tassarr,kass,d,H,Hip1,Him1,Hjp1,Hjm1);       
        end
        if(gsmvn<gsmval)
            v=vn;
            vnisgood=1;
            gsmval=gsmvn
        end
    end
    if(vnisgood==0)
        break;
    end
    
    stepsize=norm(vstep);

    end
    errvmse=sqrt(mean((v-uTarr((ittot-1)*kass+1,1:d23)).^2))
    
    itJ=modd(itJ+1,mJ);
    lastJ{itJ}=Jarr;
    clear Jarr;
    if(mJ>0)
    lastJinv{itJ}=Jinvarr;
    clear Jinvarr;
    end
    
    
    for(l=1:kass)
    l
    mnew=zeros(1,d23);
    for(p=1:Nparticles)
    vnewsim=onestep_mex(varr(p,1:d23), Tarr(1:2), 1,imax1,d,H,Hip1,Him1,Hjp1,Hjm1,BinomialArray,FactorialArray);
    vnew(p,1:d23)=vnewsim(2,1:d23);
    mnew=mnew+vnewsim(2,1:d23);
    end

    mnew=mnew/Nparticles;
    
    %inflation
    vnew=(vnew)*(1+inflation)-ones(Nparticles,1)*mnew*inflation;   
    
    %Localized way
    C=sparse(d23,d23);
    
    Vn=vnew-ones(Nparticles,1)*mnew;
    
    [CI,CJ,CV]=localizedcov2_mex(d,Nparticles,neighborsMarr,Vn,ln,neighborstrunc);

    C=sparse(CI,CJ,CV)/(Nparticles-1);

    S=spobsuvprec*C*spobsuvprec+Gammamx;
    
    Yp=Y(l+1,1:d23)'*ones(1,Nparticles)+spobsuvprec*(sigmaz*randn(d23,Nparticles)-vnew(1:Nparticles,1:d23)');
    
    varrndiff=C*spobsuvprec*(S\Yp);
    varr=vnew+varrndiff';
    
    mest=sum(varr(1:Nparticles,1:d23),1)/Nparticles;
 
    end
    C=C/(sigmaz^2);

    end
    
    if(mod(k,kass)~=0)
        ittot=nbtot+1;
        kleft=mod(k,kass);
        vfi=pushforward(imaxsim, nbblocksass, v, Tass,d,H,Hip1,Him1,Hjp1,Hjm1);
        vTarr=runshallowwater(imaxsim,nbblocksass, vfi, Tarr(1:kleft+1), kleft,d,H,Hip1,Him1,Hjp1,Hjm1);
        for(l=1:kleft+1)
            diffuv((ittot-1)*kass+l)=sqrt(mean((abs(vTarr(l,1:2*d2)-uTarr((ittot-1)*kass+l,1:2*d2))).^2));
            diffh((ittot-1)*kass+l)=sqrt(mean((abs(vTarr(l,2*d2+1:d23)-uTarr((ittot-1)*kass+l,2*d2+1:d23))).^2));
            reldiffuv((ittot-1)*kass+l)=diffuv((ittot-1)*kass+l)/sqrt(mean((uTarr((ittot-1)*kass+l,1:2*d2)).^2));
            reldiffh((ittot-1)*kass+l)=diffh((ittot-1)*kass+l)/sqrt(mean((uTarr((ittot-1)*kass+l,2*d2+1:d23)).^2));
        end              
    end
    
    ret=struct;
    ret.func='shallowwater_4dvarnew';
    ret.diffuv=diffuv;
    ret.diffh=diffh;
    ret.reldiffuv=reldiffuv;
    ret.reldiffh=reldiffh;
    ret.err0arr=err0arr;
    ret.d=d;
    ret.sigmaz=sigmaz;
    ret.Delta=Delta;
    ret.Tass=Tass;
    ret.nbblocksass=nbblocksass;
    ret.kass=kass;
    ret.nbtot=nbtot;
    ret.obsfrequv=obsfrequv;
    ret.obsfreqh=obsfreqh;
    ret.imaxJ=imaxJ;
    ret.imaxsim=imaxsim;
    ret.mJ=mJ;
    ret.prioruv=prioruv;
    ret.priorh=priorh;
    ret.pcgsteps=pcgsteps;
    ret.loc=loc;
    ret.neighborstrunc=neighborstrunc;
    ret.params=params;
    ret.runtime=toc();
end

function ret=shallowwater_4dvarnew_extend(filename,params)
kasstrue=params.kass;
kass=kasstrue;
nbblocksasstrue=params.nbblocksass;
nbblocksass=nbblocksasstrue;
mJ=params.mJ;
interpolin=params.interpolin;
neighborstrunc=params.neighborstrunc;
pcgsteps=params.pcgsteps;

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
global priormean;

global Blist4;
global Blist2;
global Blist1;
global Amx;
global interpol;
global neighborsM;
global C;


prioruv=params.prioruv;
priorh=params.priorh;


priorprecen4dvar(d);
priorprecmx0=priorprecmx;


%%%%%%%Initial definitions 
neighborsM=neighbors(neighborstrunc,d);
imaxsim = 5;
imaxJ=4;
mxextend=params.mxextend;

interpol=interpolin;
mJcall=mJ;

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
    pcgarr=zeros(1,100);
    pcgarrit=0;
    err0arr=zeros(nbtot,1);
    
    for(ittot=1:nbtot)
    ittot
    if(ittot<=mJ)nbprevJstored=ittot-1;
    else
        nbprevJstored=mJ;
    end
    
    if(ittot==1)
        [priormean,priorprecmx]=priorupdateinitialobs(priorprecmx0,priormean,Y(1,1:d23));
        v=priormean;
        %v=u0;
        err0arr(1)=sqrt(mean((v(1:2*d2)-u0(1:2*d2)).^2));
        vTarr=runshallowwater(imaxsim,nbblocksass, v, Tassarr,  kass,d,H,Hip1,Him1,Hjp1,Hjm1);
        for(l=1:kass+1)
            diffuv(l)=sqrt(mean((abs(vTarr(l,1:2*d2)-uTarr(l,1:2*d2))).^2));
            diffh(l)=sqrt(mean((abs(vTarr(l,2*d2+1:d23)-uTarr(l,2*d2+1:d23))).^2));
            reldiffuv(l)=diffuv(l)/sqrt(mean((uTarr(l,1:2*d2)).^2));
            reldiffh(l)=diffh(l)/sqrt(mean((uTarr(l,2*d2+1:d23)).^2));
        end



    else
        vfi=pushforward(imaxsim, nbblocksasstrue, v, Tass,d,H,Hip1,Him1,Hjp1,Hjm1);
        err0arr(ittot)=sqrt(mean((v(1:2*d2)-uTarr((ittot-1)*kass+1,1:2*d2)).^2));
        vTarr=runshallowwater(imaxsim,nbblocksasstrue, vfi, Tassarr, kass,d,H,Hip1,Him1,Hjp1,Hjm1);
        for(l=1:kass+1)
            diffuv((ittot-1)*kass+l)=sqrt(mean((abs(vTarr(l,1:2*d2)-uTarr((ittot-1)*kass+l,1:2*d2))).^2));
            diffh((ittot-1)*kass+l)=sqrt(mean((abs(vTarr(l,2*d2+1:d23)-uTarr((ittot-1)*kass+l,2*d2+1:d23))).^2));
            reldiffuv((ittot-1)*kass+l)=diffuv((ittot-1)*kass+l)/sqrt(mean((uTarr((ittot-1)*kass+l,1:2*d2)).^2));
            reldiffh((ittot-1)*kass+l)=diffh((ittot-1)*kass+l)/sqrt(mean((uTarr((ittot-1)*kass+l,2*d2+1:d23)).^2));
        end
        v=vfi;
        priormean=vfi;
    end

    gsmvaltrue=gsmwithprior(imaxsim,nbblocksass,uTarr((ittot-1)*kass+1,1:d23),Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Tassarr(1:kass),kass,d,H,Hip1,Him1,Hjp1,Hjm1)

    
    err0=sqrt(mean((v-uTarr((ittot-1)*kass+1,:)).^2))

    for(extendit=1:mxextend)
        stepsize=1;
        Newtonit=0;
        kass=kasstrue/mxextend*extendit;           
        nbblocksass=nbblocksasstrue/mxextend*extendit;
    while stepsize>1e-3 && Newtonit<5 && (Newtonit<1|| ittot==1)
    %tic
    
    Newtonit=Newtonit+1
    
    if(mJ==0)
        [vTass,Jarr]=runJ(imaxsim,imaxJ, v, Tassarr(1:kass), kass,d,H,Hip1,Him1,Hjp1,Hjm1);
    else
        [vTass,Jarr,Jinvarr]=runJ(imaxsim,imaxJ, v, Tassarr(1:kass), kass,d,H,Hip1,Him1,Hjp1,Hjm1);
    end
    if(Newtonit==1)
        gsmval=gsmwithprior(imaxsim,nbblocksass,v,Y(((ittot-1)*kasstrue+1):((ittot-1)*kasstrue+kass+1),1:d23),Tassarr(1:kass),kass,d,H,Hip1,Him1,Hjp1,Hjm1);
    end
    %gsmJ(v,vTass,Y((ittot-1)*kass+1:ittot*kass+1,1:d23))

    gradnJ=gradJ(v,vTass,Y(((ittot-1)*kasstrue+1):((ittot-1)*kasstrue+kass+1),1:d23),Jarr,kass,d);
      
  
    
    kcall=kass;dcall=d;Jarrcall=Jarr;
    %vstep=(pcg(@HessJcallnew,gradnJ,1e-4,pcgsteps))';
    if(ittot==1)
        vstep=(pcg(@HessJcall,gradnJ,1e-6,pcgsteps))';
        %[vstep,flag,rr,pcgiter] =(pcg(@HessJcall,gradnJ,1e-6,pcgsteps));
    else
    %[vstep,flag,rr,pcgiter] =(pcg(@HessJcallnew,gradnJ,1e-6,pcgsteps));
    vstep=(pcg(@HessJcallnew,gradnJ,1e-6,pcgsteps))';
    end
    


    vnisgood=0;
    vstep=2*vstep;
    vnit=0;
    while(vnisgood==0 && vnit<=10)
        vnit=vnit+1;
        vstep=vstep/2;
        vn=v-vstep;
        gsmvn=gsmwithprior(imaxsim,nbblocksass,vn,Y((ittot-1)*kass+1:ittot*kass+1,1:d23),Tassarr,kass,d,H,Hip1,Him1,Hjp1,Hjm1);
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
    end
    
    errvmse=sqrt(mean((v-uTarr((ittot-1)*kass+1,1:d23)).^2))
    
    itJ=modd(itJ+1,mJ);
    lastJ{itJ}=Jarr;
    clear Jarr;
    if(mJ>0)
    lastJinv{itJ}=Jinvarr;
    clear Jinvarr;
    end

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
    ret.func='shallowwater_4dvarnew_extend';
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
    ret.neighborstrunc=neighborstrunc;
    ret.pcgarr=pcgarr;
    ret.params=params;
    ret.runtime=toc();
end

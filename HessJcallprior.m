function Hessw=HessJcallprior(w)
global kcall;
global dcall;
global Jarrcall;

global priorprecmx;
global priormean;
global spobsuvprec;

global lastJ;
global lastJinv;
global itJ;
global nbprevJstored;
global mJcall;
global interpol;


if(nbprevJstored==0)
    Hessw=priorprecmx*w;return;
end

mind=zeros(nbprevJstored,1);
if(nbprevJstored<mJcall)
    mind=nbprevJstored:-1:1;
else
for(it=1:mJcall)
    mind(it)=modd(itJ-it+1,mJcall);
end
end

d23=dcall^2*3;


warr=zeros(d23,nbprevJstored);

warr(1:d23,1)=w;

for(it=kcall:-1:1)
%    warr(1:d23,1)=lastJinv{mind(1)}{it}*warr(1:d23,1);
if(it==kcall|| it==1|| mod(it+1,interpol)==0) warr(1:d23,1)=lastJinv{mind(1)}{it}*warr(1:d23,1);
else
    itless=max(1,it-mod(it,interpol)-1);
    itmore=min(kcall,it+interpol-mod(it,interpol)-1);
    lambda=(it-itless)/(itmore-itless);
    warr(1:d23,1)=lastJinv{mind(1)}{itless}*warr(1:d23,1)*(1-lambda)+lastJinv{mind(1)}{itmore}*warr(1:d23,1)*lambda;
end

end
%nw=norm(warr(1:d23,1))


for(itJac=2:nbprevJstored)
    warr(1:d23,itJac)=warr(1:d23,itJac-1);
for(it=kcall:-1:1)
%    warr(1:d23,itJac)=lastJinv{mind(itJac)}{it}*warr(1:d23,itJac);
if(it==kcall|| it==1|| mod(it+1,interpol)==0) warr(1:d23,itJac)=lastJinv{mind(itJac)}{it}*warr(1:d23,itJac);
else
    itless=max(1,it-mod(it,interpol)-1);
    itmore=min(kcall,it+interpol-mod(it,interpol)-1);
    lambda=(it-itless)/(itmore-itless);
    warr(1:d23,itJac)=lastJinv{mind(itJac)}{itless}*warr(1:d23,itJac)*(1-lambda)+lastJinv{mind(itJac)}{itmore}*warr(1:d23,itJac)*(lambda);
end
end
end

Hessw=priorprecmx*warr(1:d23,nbprevJstored)+HessJnoprior(warr(1:d23,nbprevJstored),lastJ{mind(nbprevJstored)},kcall,dcall);



for(it=1:kcall)
    %Hessw=lastJinv{mind(nbprevJstored)}{it}'*Hessw;
    if(it==kcall|| it==1|| mod(it+1,interpol)==0) Hessw=lastJinv{mind(nbprevJstored)}{it}'*Hessw;
    else
        itless=max(1,it-mod(it,interpol)-1);
        itmore=min(kcall,it+interpol-mod(it,interpol)-1);
        lambda=(it-itless)/(itmore-itless);
        Hessw=lastJinv{mind(nbprevJstored)}{itless}'*Hessw*(1-lambda)+lastJinv{mind(nbprevJstored)}{itmore}'*Hessw*(lambda);
    end        

end

for(itJac=nbprevJstored-1:-1:1)
    Hessw=Hessw+HessJnoprior(warr(1:d23,itJac),lastJ{mind(itJac)},kcall,dcall);
    for(it=1:kcall)
        %Hessw=lastJinv{mind(itJac)}{it}'*Hessw;
        if(it==kcall|| it==1|| mod(it+1,interpol)==0) Hessw=lastJinv{mind(itJac)}{it}'*Hessw;
        else
            itless=max(1,it-mod(it,interpol)-1);
            itmore=min(kcall,it+interpol-mod(it,interpol)-1);
            lambda=(it-itless)/(itmore-itless);
            Hessw=lastJinv{mind(itJac)}{itless}'*Hessw*(1-lambda)+lastJinv{mind(itJac)}{itmore}'*Hessw*lambda;
        end        
    end
end
%nHw=norm(Hessw)

end
function covinteract2=pmx_plot(ret)
pmx=ret.pmx;
cov=inv(pmx);
d=ret.d;
d23=d^2*3;
covinteract=zeros(d+1,1);
covit=zeros(d+1,1);
for(i1=1:d)
    for(j1=1:d)
        for(i2=1:d)
            for(j2=1:d)
                di=min(abs(i1-i2),d-abs(i1-i2));
                dj=min(abs(j1-j2),d-abs(j1-j2));
                dij=di+dj;
                covit(dij+1)=covit(dij+1)+1;

                for(uvh1=1:3)
                    for(uvh2=1:3)               
                        ij1=(i1-1)*d+j1+(uvh1-1)*d^2;
                        ij2=(i2-1)*d+j2+(uvh2-1)*d^2;
                        covinteract(dij+1)=covinteract(dij+1)+abs(cov(ij1,ij2));
                    end
                end
            end
        end
    end
end
covinteract2=covinteract./covit/9;
covinteract2=covinteract2/covinteract2(1);
end
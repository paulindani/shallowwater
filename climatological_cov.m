function [Cmx,Cmean]=climatological_cov(inputfile)
load(inputfile)

Cmx=zeros(d23,d23);
for(i=1:d23)
    for(j=1:i)      
    Cmx(i,j)=mean(uTarr(:,i).*uTarr(:,j));
    Cmx(j,i)=Cmx(i,j);
    end
end
%This step adds a small identity matrix in order to avoid invertibility
%issues

Cmx=Cmx+norm(Cmx)/100 * eye(d23);

Cmean=mean(uTarr,1);
end
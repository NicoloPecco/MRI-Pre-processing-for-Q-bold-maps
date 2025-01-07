function [amap,cmap,devmap]=imagefit2param(zeit,serie,x,y,g1,g2,p,noisefakt,maxdev)
% function [amap,cmap]=imagefit2param(zeit,serie,x,y,g1,g2,p,noisefakt,maxdev)
% 2 parameter fit:  wert[i] = a*{x+y*exp(-zeit[i]/c)}
% zeit: vector containing time points
% serie: 3D set containing weighted images
% x,y according to fit formula
% g1,g2: lower/upper limit for c
% p: accuracy in relative values
% noisefakt: noise factor
% maxdev: maximal deviation of original from fitted values, in percent


[nphase,nread,nslc]=size(serie);
zeit=reshape(zeit,1,nslc);
%serie=abs(serie);

noise=serie(1:nphase/8,1:nread/8,1);
noise=reshape(noise,1,nphase/8*nread/8);
noise=std(noise(:))/0.655;

maske=zeros(nphase,nread);
amap=zeros(nphase,nread);
cmap=zeros(nphase,nread);
devmap=zeros(nphase,nread);

for i=1:nslc
 maske=maske+(serie(:,:,i)>noisefakt*noise);

end

for i=1:nphase
    for k=1:nread
        if maske(i,k)>0
            wert=squeeze(serie(i,k,:));
            wert=reshape(wert,1,nslc);
            [a,c,dev]=fit2param(zeit,wert,x,y,g1,g2,p); 
            if (dev<maxdev && c<g2*(1-2*p))
                amap(i,k)=a;
                cmap(i,k)=c;
                devmap(i,k)=dev;
            end
        end
    end
end

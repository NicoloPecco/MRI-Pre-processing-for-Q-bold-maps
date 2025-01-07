function [a,c,dev]=fit2param(zeit,wert,x,y,g1,g2,p)
% function [a,c,dev]=fit2param(zeit,wert,x,y,g1,g2,p
% 2 parameter fit:  wert[i] = a*{x+y*exp(-zeit[i]/c)}
% zeit: vector containing time points
% wert: vector containing values
% x,y according to fit formula
% g1,g2: lower/upper limit for c
% p: accuracy in relative values
% dev: deviation of original from fitted values, in percent

num=length(zeit);

fahne=1;
while fahne==1
 c=(g1+g2)/2;
 hilf=exp(-zeit./c);
 s1=sum(hilf);
 s2=sum(hilf.*hilf);
 s3=sum(wert.*hilf);
 s4=sum(wert);
 s5=sum(zeit.*hilf);
 s6=sum(zeit.*hilf.*hilf);
 s7=sum(zeit.*wert.*hilf);
 a=(x*s4+y*s3)/(num*x*x+2*x*y*s1+y*y*s2);
 slopp=2*a*y/(c*c)*(s5*a*x+s6*a*y-s7);
 if (g2-g1)/c < p
  fahne=0;
 end
 if slopp>=0
  g2=c;
 end
 if slopp<0
  g1=c;
 end
end

hilf=a*x+a*y*exp(-zeit./c);
s8=sum(abs(wert-hilf));
s9=sum(abs(wert));
dev=100*s8/s9;

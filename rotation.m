% Script written by:
% Zhuo Li (zhuol7@student.unimelb.edu.au)
% The University of Melbourne

function [xf,yf,zf]=rotation(xi,yi,zi,R)

 
I=size(xi,1);
J=size(xi,2);
 
xf=zeros(I,J);
yf=zeros(I,J);
zf=zeros(I,J);
 
for ii=1:I
	for jj=1:J
    	vector=[xi(ii,jj);yi(ii,jj);zi(ii,jj)];
    	vector=R*vector;
        	xf(ii,jj)=vector(1);
        	yf(ii,jj)=vector(2);
        	zf(ii,jj)=vector(3);
	end
end

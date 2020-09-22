function Gout=Leffect(f,z1,z2,wavelength,dx,N,I)
%Simulation of Fourier transformation of bitmap images
% clear
%%  read bmp figure  preparation
M=-1;
% f=100; M=-6; %um 
% z1=(1-1/M)*f;
% z2=(1-M)*f; 
% z1=f;z2=f;
% wavelength=0.633; %微米
k=2*pi/wavelength;
% dx=0.3; N=300;%input gap%300*300像素，每个像素的大小为300nm*300nm
x0=((1:N)-0.5*N)*dx; y0=x0;
[xx0,yy0]=meshgrid(x0,y0);
r0=sqrt(xx0.^2+yy0.^2);
fai=atan2(yy0,xx0);
% pupil=zeros(size(r));%define P
% pupil(r<=5)=1;%r=10μm
% figure(1);imagesc(pupil);title('circle');
% I=imread('a.bmp','bmp'); %Input bitmap image 
% I=I(:,:,1);
figure(1);imagesc(x0,y0,I);title('input star'); %displaying input
colormap(gray(255));
% image(I)
% axis off
%% first R-S diffraction

x=x0;y=y0;%x,y of lens 
[xx,yy]=meshgrid(x,y);
r=sqrt(xx.^2+yy.^2);
fai=atan2(yy,xx);

Ei=I;
Eout1=zeros(N,N);
Eout1(:,:)=RS2D(x0,y0,Ei,x,y,wavelength,z1);
figure(2);imagesc(abs(Eout1));title('middle layer'); %displaying input
%% lens effect

Eoutp=Eout1;
Eoutp(r>=420)=0;%define P of lens %r=30μm

lens=exp(-1i*k/2/f*(r.^2));

I2=Eoutp.*lens;%new output

%% second R-S diffraction

dxi=0.3; Ni=-M*N;%input gap%300*300像素，每个像素的大小为300nm*300nm
xi=((1:Ni)-0.5*Ni)*dxi; yi=xi;
[xxi,yyi]=meshgrid(xi,yi);
% xi=-M*x0;yi=-M*y0;%x,y of lens 
% [xxi,yyi]=meshgrid(xi,yi);
r=sqrt(xxi.^2+yyi.^2);
faii=atan2(yyi,xxi);

% NR=6*N;
% xreference=linspace(-300,300,NR);
% yreference=linspace(-300,300,NR);

Ei=I2;
Gout=zeros(Ni,Ni);
Gout(:,:)=abs(RS2D(x,y,Ei,xi,yi,wavelength,z2));
figure(3);imagesc(xi,yi,abs(Gout));title('output'); %displaying input
end


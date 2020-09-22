%% ????????????
clear
%%  read bmp figure  preparation
Rmin=10;Rmax=60;
HRmin=10;HRmax=40;
f=100; M=-1; %um 
% z1=(1-1/M)*f;
% z2=(1-M)*f; 
z1=f;z2=2*f;
wavelength=0.633; %微米
k=2*pi/wavelength;
dx=0.3; N=300;%input gap%300*300像素，每个像素的大小为300nm*300nm
x0=((1:N)-0.5*N)*dx; y0=x0;
[xx0,yy0]=meshgrid(x0,y0);
r0=sqrt(xx0.^2+yy0.^2);
fai=atan2(yy0,xx0);
% pupil=zeros(size(r));%define P
% pupil(r<=5)=1;%r=10μm
% figure(1);imagesc(pupil);title('circle');
%I=imread('a','bmp'); %Input bitmap image 
I = imread('a','bmp'); 
I=imresize(I,[300 300]);
I=I(:,:,1);
figure(1);imagesc(x0,y0,I);title('input image'); %displaying input
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
% figure(2);imagesc(abs(Eout1));title('middle layer'); %displaying input
%% 1st lens effect

Eoutp1=Eout1;
% Eoutp1(r<=Rmin)=0;%define P of lens %r=30μm
Eoutp1(r>=Rmax)=0;

lens=exp(-1i*k/2/f*(r.^2));

I2=Eoutp1.*lens;%new output

%% second R-S diffraction

dxi=0.3; Ni=-M*N;%input gap%300*300像素，每个像素的大小为300nm*300nm
xi=((1:Ni)-0.5*Ni)*dxi; yi=xi;
[xxi,yyi]=meshgrid(xi,yi);
% xi=-M*x0;yi=-M*y0;%x,y of lens 
% [xxi,yyi]=meshgrid(xi,yi);
r=sqrt(xxi.^2+yyi.^2);
faii=atan2(yyi,xxi);

Ei=I2;
Eouti2=zeros(Ni,Ni);
Eouti2(:,:)=RS2DH(x,y,Ei,xi,yi,wavelength,z1);
% figure(3);imagesc(xi,yi,abs(Eouti2));title('second R-S diffraction output'); %displaying input

%% H's effect
H=r.*0+1;
%H(r<=HRmin)=0; %????????????+?????
%H(r>=HRmax)=0;
%H=-4*pi^2*((xx0.*(wavelength*f)^(-1))^2+(yy0.*(wavelength*f)^(-1))^2);%laplace
H=-1i*2*pi*((xx0.*(wavelength*f)^(-1))^1+(yy0.*(wavelength*f)^(-1))^1);%only one derivation
%x=zeros(1,150);
%y=zeros(1,150);
%H(150,150)=0;H(151,150)=0;H(150,151)=0;H(149,150)=0;H(150,149)=0;H(151,151)=0;H(151,149)=0;H(149,151)=0;H(149,149)=0;
%H(152,150)=1;H(151,152)=1;H(152,151)=1;H(148,150)=1;H(150,148)=1;H(152,152)=1;H(152,148)=1;H(148,152)=1;H(148,148)=1;
% lens=exp(-1i*k/2/f*(r.^2));
EoutH=Eouti2.*H;
%EoutH=Eouti2;
% I3=Eoutp2.*lens;%new output

%% third R-S diffraction

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

Ei=EoutH;
Eouti3=zeros(Ni,Ni);
Eouti3(:,:)=RS2DH(x,y,Ei,xi,yi,wavelength,z1);
% figure(3);imagesc(xi,yi,abs(Eouti3));title('second R-S diffraction output'); %displaying input

% FI=fft2(I);
% FI=fftshift(FI);
% max1=max(FI);
% max2=max(max1);
% scale=1.0/max2;
% FI=FI.*scale;
% figure(2) %Gray scale image of the absolute value of transform

%% 2nd lens effect

Eoutp2=Eouti3;
% Eoutp2(r<=Rmin)=0;%define P of lens %r=30μm
Eoutp2(r>=Rmax)=0;

lens=exp(-1i*k/2/f*(r.^2));

I3=Eoutp2.*lens;%new output

%% forth R-S diffraction

% dxi=0.3; Ni=-M*N;%input gap%300*300像素，每个像素的大小为300nm*300nm
% xi=((1:Ni)-0.5*Ni)*dxi; yi=xi;
% [xxi,yyi]=meshgrid(xi,yi);
% xi=-M*x0;yi=-M*y0;%x,y of lens 
% [xxi,yyi]=meshgrid(xi,yi);
% r=sqrt(xxi.^2+yyi.^2);
% faii=atan2(yyi,xxi);

% NR=6*N;
% xreference=linspace(-300,300,NR);
% yreference=linspace(-300,300,NR);

Ei=I3;
Eouti4=zeros(Ni,Ni);
Eouti4(:,:)=RS2D(x,y,Ei,xi,yi,wavelength,z1);
figure(3);imagesc(xi,yi,abs(Eouti4));title('output from 4F system'); %displaying input

%fft2Dbitmap_image.m
%Simulation of Fourier transformation of bitmap images
clear
%%  read bmp figure  preparation

f=100; M=-3; %um 
z1=(1-1/M)*f;
z2=(1-M)*f; 
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
I=imread('a.bmp','bmp'); %Input bitmap image 
I=I(:,:,1);
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
xi=linspace(-300,300,Ni);yi=xi;
% yi=linspace(-300,300,Ni);
% xi=((1:Ni)-0.5*Ni)*dxi; yi=xi;
[xxi,yyi]=meshgrid(xi,yi);
% xi=-M*x0;yi=-M*y0;%x,y of lens 
% [xxi,yyi]=meshgrid(xi,yi);
r=sqrt(xxi.^2+yyi.^2);
faii=atan2(yyi,xxi);

% NR=6*N;
% xreference=linspace(-300,300,NR);
% yreference=linspace(-300,300,NR);

Ei=I2;
Eouti=zeros(Ni,Ni);
Eouti(:,:)=RS2D(x,y,Ei,xi,yi,wavelength,z2);
figure(3);imagesc(xi,yi,abs(Eouti));title('output'); %displaying input

% FI=fft2(I);
% FI=fftshift(FI);
% max1=max(FI);
% max2=max(max1);
% scale=1.0/max2;
% FI=FI.*scale;
% figure(2) %Gray scale image of the absolute value of transform


function Eout = RS2D(x,y,Ei,xo,yo,wavelength,z) %xo,yo=output
k=2*pi/wavelength;
N=length(x);
M=length(xo);
dx=abs(x(2)-x(1));
X=zeros(N+M-1,1);%extend x->X
Y=X;
X(1:N-1)=xo(1)-x(N:-1:2);
X(N:N+M-1)=xo(1:M)-x(1);
Y(1:N-1)=yo(1)-y(N:-1:2);
Y(N:N+M-1)=yo(1:M)-y(1);
[XX YY]=meshgrid(X,Y);
r=sqrt(XX.^2+YY.^2+z^2);

H=1/2/pi.*exp(1i*k*r)./r.^2*z.*(1./r-1i*k);
U=zeros(M+N-1);
U((1:N)-round(N/2)+round((N+M-1)/2),(1:N)-round(N/2)+round((N+M-1)/2))=Ei;%Ei at the centre of U

e1=(fft2(U)); 
e2=(fft2(H));%FFT 2demensions
e=fftshift(ifft2((e1.*e2)))*dx^2;%fftshift modify fft2's side effect
Eout=e(-1+(1:M)-1*round(M/2)+1*(round((M+N-1)/2)),-1+(1:M)-1*round(M/2)+1*round((N+M-1)/2));
end



%% target
A = imread('target.jpg','jpg'); %Input bitmap image 
lengthr=96/25.4; hightr=96/25.4; %para of fig
B = imresize(A,[500*hightr  500*lengthr]); %实际500mm,此处认为500μm
Target = rgb2gray(B);
imagesc(Target);
[TM,TN]=size(Target);%size of the target
%% parameter
lambda=632.8e-3;
radius=100;% μm
Z = 400;%distance
% D = 500;%size of the target
dx = 0.3;  %length of each pixel
x1= -((TM-1)*dx/2):dx:((TM-1)*dx/2); y1 = -((TN-1)*dx/2):dx:((TN-1)*dx/2); % result's axis

%% inputs
x0=x1; y0=y1; %inputs' geo size=500*500
[xx0,yy0]=meshgrid(x0,y0);
%source
f=1e3;% source's para?
Source = exp(-1i/2*(xx0.^2+yy0.^2)/f^2); 
r0=sqrt(xx0.^2+yy0.^2);
Source(r0>=radius)=0;
imagesc(abs(Source))%source's image
%phase
%Phi = 0*(xx0.^2+yy0.^2)+1+1i;
%Phi(r0>=radius)=0;
%B = abs(Source) .* exp(1i*angle(Phi));%transmission input

%% transmission 
% A = fftshift(ifft2(fftshift(Target)));
epsilon=1e8;%convergent index
RMSE=1e10;index=0;
Target=double(Target);
Phi = fftshift(ifft2(fftshift(Target)));%phase
while RMSE>epsilon && index<10
    B = abs(Source) .* exp(1i*angle(Phi));%transmission input
%     Result = fftshift(fft2(fftshift(B)));
    Result = RS2D(x0,y0,B,x1,y1,lambda,Z);%result 
    D = abs(Target) .* exp(1i*angle(Result));
    %D = abs(Target) .* cos(angle(Result)*180/pi)+1i*abs(Target) .* sind(angle(Result)*180/pi); % optimize by combining angle of C and amplification of Target
%     Phi = fftshift(ifft2(fftshift(D))); %IFFT
    Phi=iRS(x1,y1,D,x0,y0,lambda,Z);
    
    imagesc(abs(Result)) %Present current pattern
    index=index+1;
    title(sprintf('%d',index));
    RMSE=sum((sum(abs(Result-Target)))');
    pause(0.5)
end


%% RayleighSommerfeld
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
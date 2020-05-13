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
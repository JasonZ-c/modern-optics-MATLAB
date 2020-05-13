[X,map]=imread('a.bmp');
Y=zeros(300,300);
%each pixel=300nm*300nm
x0=-44.85:0.3:44.85;
y0=-44.85:0.3:44.85;
x0=1e-6*x0;
y0=1e-6*y0;

%object's info
U0=zeros(length(x0),length(y0));
for i=1:length(x0)
    for j=1:length(y0)
        if X(i,j)>10    
            U0(i,j)=1;
        else
            U0(i,j)=0;
        end
    end
end
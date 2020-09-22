close all;  
clear;
clc;  
LCP=load('Tansmission of LCP.mat');
RCP=load('Tansmission of RCP.mat');
T_all=load('Tansmission all.mat');

[mL,pL]=max(max(max(LCP(1).TransLCP)));
a_max=max(max(max(LCP(1).TransLCP)));
[x,y,z]=find(LCP(1).TransLCP==max(max(max(LCP(1).TransLCP))));
[mR,pR]=max(max(max(RCP(1).TransRCP)));
[xR,yR,zR]=find(RCP(1).TransRCP==max(max(max(RCP(1).TransRCP))));

[mA,pA]=max(max(max(T_all(1).TransT)));
[xA,yA,zA]=find(RCP(1).TransRCP==max(max(max(T_all(1).TransT))));
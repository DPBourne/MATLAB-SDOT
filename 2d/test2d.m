x0=0;y0=0;
l1=1;l2=2;

Omega=[x0,y0,x0+l1,y0+l2];

N=500000;
X=[x0 y0]+rand(N,2)*diag([l1,l2]);
w0=0.01*rand(N,1)/N;

periodic=true;

percent_tol=0.1;
m=ones(N,1)*l1*l2/N;

tic
[w,percent_err,v,EXITFLAG] = SDOT2d_damped_Newton(w0,X,m,Omega,periodic,percent_tol);
toc

% $$$ [g0,dg0,d2g0]=kantorovich2d(w0,X,m,Omega,periodic);
% $$$ 
% $$$ for j=1:N,
% $$$     wp=w0;wp(j)=wp(j)+epsilon;
% $$$     wm=w0;wm(j)=wm(j)-epsilon;
% $$$ 
% $$$     [gp,dgp,d2gp]=kantorovich2d(wp,X,m,Omega,periodic);
% $$$     [gm,dgm,d2gm]=kantorovich2d(wm,X,m,Omega,periodic);
% $$$ 
% $$$     approx_d2g(:,j)=(dgp-dgm)/2/epsilon;
% $$$ end
% $$$ 
% $$$ diff=d2g0-approx_d2g;


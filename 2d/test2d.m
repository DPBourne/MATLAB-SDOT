x0=0;y0=0;
l1=1;l2=2;

Omega=[x0,y0,x0+l1,y0+l2];

N=5000;
X=[x0 y0]+rand(N,2)*diag([l1,l2]);
w0=zeros(N,1);

periodic=true;

percent_tol=0.1;
m=ones(N,1)*l1*l2/N;

tic
[w,percent_err,v,EXITFLAG] = SDOT2d_damped_Newton(w0,X,m,Omega,periodic,percent_tol);
toc

function report = run_benchmarks()
%RUN_BENCHMARKS Analytic diffusion, advection and local reaction benchmarks.
% Uses the production five-field RHS and stage integrator in every run.
root=fileparts(fileparts(mfilename('fullpath')));report=struct('runtime',version);
ns=[10 20 40 80 160];report.diffusion=zeros(numel(ns),5);
for j=1:numel(ns)
 p=base(ns(j));p.L=.1;p.ks=.2;p.rtol=1e-10;p.atol=1e-12;
 C=p.rho_bd*(p.cd+p.cl*p.X0);alpha=p.ks/C;tf=.2*p.L^2/alpha;
 z=((1:p.N)'-.5)*p.L/p.N;am=5*sin(pi/(2*p.N))/(pi/(2*p.N));
 y=plt.initial(p);y(p.N+1:2*p.N)=20+am*cos(pi*z/p.L);
 o=struct('tEnd',tf,'initialState',y,'maxStep',tf/30,'target',60);
 out=plt.simulate(p,plt.schedule(0,[20 0 0]),o);
 exact=20+am*cos(pi*z/p.L)*exp(-alpha*(pi/p.L)^2*tf);
 err=out.y(end,p.N+1:2*p.N)'-exact;
 report.diffusion(j,:)=[p.N,sqrt(mean(err.^2)),max(abs(err)),out.elapsedSeconds,tf];
 fprintf('diffusion N=%d L2=%.5g\n',p.N,report.diffusion(j,2));
end
report.advection=zeros(numel(ns),5);
for j=1:numel(ns)
 p=base(ns(j));p.Yin=.1;p.rtol=1e-10;p.atol=1e-12;u=.01;
 eg=p.eps-p.rho_bd*p.X0/p.rho_l;v=u/eg;tf=.15*p.L/v;
 faces=(0:p.N)'*p.L/p.N;z=.5*(faces(1:end-1)+faces(2:end));
 [~,rd]=plt.inlet([20 0 u],p);Y=p.Yin+.05*cellBump(faces,0);
 y=plt.initial(p);y(4*p.N+1:5*p.N)=eg*rd*Y;
 o=struct('tEnd',tf,'initialState',y,'maxStep',tf/60,'target',60);
 out=plt.simulate(p,plt.schedule(0,[20 0 u]),o);
 exact=p.Yin+.05*cellBump(faces,v*tf);
 pred=out.y(end,4*p.N+1:5*p.N)'/(eg*rd);err=pred-exact;
 report.advection(j,:)=[p.N,mean(abs(err)),max(abs(err)),out.elapsedSeconds,tf];
 if j==numel(ns),report.advectionProfile=[z,exact,pred];end
 fprintf('advection N=%d L1=%.5g\n',p.N,report.advection(j,2));
end
rtols=[1e-3 1e-5 1e-7 1e-9];report.reaction=zeros(numel(rtols),4);
for j=1:numel(rtols)
 p=base(4);p.bio=true;p.Ho=0;p.rtol=rtols(j);p.atol=rtols(j)*1e-3;
 eg=p.eps-p.rho_bd*p.X0/p.rho_l;[~,rd]=plt.inlet([20 0 0],p);
 p.user.k=2e-4;p.user.storage=eg*rd;p.reactionFcn=@decay;
 tf=3/p.user.k;y=plt.initial(p);
 out=plt.simulate(p,plt.schedule(0,[20 0 0]),struct('tEnd',tf,'target',60,'maxStep',tf));
 exact=y(4*p.N+1)*exp(-p.user.k*out.t);
 pred=out.y(:,4*p.N+1);err=max(abs(pred-exact))/y(4*p.N+1);
 report.reaction(j,:)=[rtols(j),err,out.elapsedSeconds,numel(out.t)];
 fprintf('reaction tol=%g error=%.5g\n',rtols(j),err);
end
report.diffusionOrder=log(report.diffusion(1:end-1,2)./report.diffusion(2:end,2))/log(2);
report.advectionOrder=log(report.advection(1:end-1,2)./report.advection(2:end,2))/log(2);
save(fullfile(root,'results','benchmark_report.mat'),'report','-v7');
assert(all(report.diffusionOrder>1.95&report.diffusionOrder<2.05));
assert(report.advectionOrder(end)>.8&&report.advectionOrder(end)<1.2);
assert(report.reaction(end,2)<1e-7);
save(fullfile(root,'results','benchmark_report.mat'),'report','-v7');
end
function p=base(n)
p=plt.defaults();p.N=n;p.T0=20;p.Tamb=20;p.RH0=0;
p.phase=false;p.bio=false;p.Uw=0;p.Dw=0;p.ks=0;p.kg=0;p.Dv=0;p.Do=0;
p.exchangeFcn=@noexchange;
end
function [h,S]=noexchange(s,~),h=zeros(size(s.Tg));S=h;end
function R=decay(~,~,Y,p),R=p.user.k*p.user.storage*Y;end
function a=cellBump(faces,shift)
% Exact cell average of sin^4 compact pulse on [0.10+shift,0.70+shift].
w=.60;x=min(w,max(0,faces-shift-.10));
I=3*x/8-w*sin(2*pi*x/w)/(4*pi)+w*sin(4*pi*x/w)/(32*pi);
a=diff(I)./diff(faces);
end

function report = test_porousltne(fullRun)
%TEST_POROUSLTNE Executable invariants, independent fixtures and trajectory tests.
% test_porousltne(false) omits the 10/20/40/60-cell grid suite.
if nargin<1,fullRun=true;end
root=fileparts(fileparts(mfilename('fullpath')));fdir=fullfile(root,'tests','fixtures');
report=struct('version','0.2.0','runtime',version,'passed',{{}},'rhsError',0);
p=plt.defaults();p.N=6;plt.validate(p);
assert(abs(plt.psat(0)-611.21)<1e-10);pass('saturation at zero Celsius');
assert(abs(plt.quantile((1:10)',.05)-1.45)<1e-14);pass('explicit quantile convention');
assert(all(plt.diffFlux([1;2;3],0,1)==0));pass('zero diffusivity');
for j=0:2
 f=load(fullfile(fdir,sprintf('rhs_%d.mat',j)));r=plt.rhs(0,f.y,f.control,p);
 e=max(abs(r-f.rhs)./max(1,abs(f.rhs)));report.rhsError=max(report.rhsError,e);
 assert(e<1e-11,'RHS fixture mismatch.');
 % Independent finite-difference check of total enthalpy derivative.
 h=1e-4;[E1,W1]=plt.inventory(f.y+h*r,p);[E2,W2]=plt.inventory(f.y-h*r,p);
 assert(abs((E1-E2)/(2*h)-sum(r(end-6:end-4)))<2e-7);
 assert(abs((W1-W2)/(2*h)-r(end-3))<2e-7);
end
pass('three Python RHS fixtures');pass('pointwise energy and water conservation');
q=p;q.bio=false;q.phase=false;q.T0=q.Tamb;q.RH0=q.RHamb;
r=plt.rhs(0,plt.initial(q),[q.T0 q.RH0 .03],q);
assert(max(abs(r(1:5*p.N)))<1e-10);pass('uniform equilibrium');
q=p;q.N=1;reject(@()plt.validate(q));
q=p;q.eps=.1;reject(@()plt.validate(q));
reject(@()plt.schedule([0 0],[30 .9 .03;30 .9 .03]));
reject(@()plt.schedule(0,[30 1.2 .03]));pass('invalid inputs rejected');
q=p;q.reactionFcn=@(T,X,Y,P) zeros(size(T));
s=plt.stateProperties(plt.initial(q),[30 .95 .03],q);assert(all(s.R==0));
q=p;q.exchangeFcn=@exchange; s=plt.stateProperties(plt.initial(q),[30 .95 .03],q);
assert(all(s.hvol==7)&&all(s.Sp==0));pass('constitutive callbacks');
o=struct('tEnd',600,'maxStep',30);
a=plt.simulate(p,plt.schedule(0,[30 .95 .03]),o);
b=plt.simulate(p,plt.schedule([0 300],[30 .95 .03;30 .95 .03]),o);
assert(max(abs(a.y(end,1:12)-b.y(end,1:12)))<2e-3);pass('stage restart continuity');
assert(isempty(a.eventTime));pass('unreached threshold');
q=p;q.T0=25;q.Tamb=25;
a=plt.simulate(q,plt.schedule(0,[25 .95 .03]),struct('stopAtTarget',true));
assert(a.eventTime==0&&numel(a.t)==1);pass('initially reached threshold');
q.reactionFcn=@(T,X,Y,P) zeros(size(T));
a=plt.simulate(q,plt.schedule(0,[25 .95 .03]),struct('stopAtTarget',true));
tmp=tempname;mkdir(tmp);cleaner=onCleanup(@()rmdir(tmp,'s'));plt.export(a,tmp);
saved=load(fullfile(tmp,'result.mat'));assert(ischar(saved.out.parameters.reactionFcn));
assert(isequal(saved.out.y,a.y));pass('portable callback export');
if fullRun
 report.grid=zeros(4,8);ns=[10 20 40 60];
 for i=1:4
  q=plt.defaults();q.N=ns(i);out=plt.simulate(q,plt.schedule(0,[30 .95 .03]));
  ref=load(fullfile(fdir,sprintf('python_grid%d.mat',q.N)));
  dt=abs(out.metrics.activationHours-ref.activationHours);
  sampleTimes=(0:4:24)'*3600;
  u=interp1(out.t,out.y(:,q.N+1:2*q.N),sampleTimes);
  v=interp1(ref.t,ref.Ts,sampleTimes);
  e=max(abs(u(:)-v(:)));
  assert(dt<0.01&&e<0.02,'Trajectory parity failed.');
  assert(out.metrics.energyResidualPercent<0.02&&out.metrics.waterResidualPercent<0.02);
  report.grid(i,:)=[q.N,out.metrics.activationHours,dt,e, ...
   out.metrics.energyResidualPercent,out.metrics.waterResidualPercent, ...
   out.metrics.conditioningKWh,out.elapsedSeconds];
  plt.export(out,fullfile(root,'results',sprintf('grid%d',q.N)));
  fprintf('N=%d t95=%.8f h parity=%.4g h temperature=%.4g C\n',q.N, ...
   out.metrics.activationHours,dt,e);
 end
 pass('four-grid integration and Python trajectory parity');
end
q=p;q.Dv=0;q.Do=0;q.kg=0;q.ks=0;q.Ho=0;plt.validate(q);
z=plt.rhs(0,plt.initial(q),[30 .95 .03],q);assert(all(isfinite(z)));
s=plt.stateProperties(plt.initial(q),[30 .95 .03],q);assert(all(s.Sp==0));
pass('zero transport coefficients');
q=p;q.user.rate=2;plt.validate(q);q.user.rate=NaN;reject(@()plt.validate(q));
pass('numeric user constants');
for label={'inert_granulate','wet_granulate','reactive_bed'},plt.material(label{1});end
pass('material presets');
f=load(fullfile(fdir,'rhs_1.mat'));q=p;z=plt.rhs(0,f.y,f.control,q);
s=plt.stateProperties(f.y,f.control,q);dz=q.L/q.N;
Fo=s.gda*[q.Yin;s.Y]+plt.diffFlux(s.Y,q.Do*s.rd.*s.eg.^(1+q.nO),dz,q.Yin);
assert(abs(sum(z(4*q.N+1:5*q.N))*dz-(Fo(1)-Fo(end)-sum(s.R)*dz))<1e-12);
pass('pointwise oxygen conservation');
f=load(fullfile(fdir,'rhs_1.mat'));J=plt.jacobian(0,f.y,f.control,p);
baseRhs=plt.rhs(0,f.y,f.control,p);sc=kron([1;1;1e-3;1;1e-2],ones(p.N,1));
for col=1:5*p.N
 h=sqrt(eps)*max(abs(f.y(col)),sc(col));yp=f.y;yp(col)=yp(col)+h;
 exact=(plt.rhs(0,yp,f.control,p)-baseRhs)/h;
 assert(max(abs(full(J(:,col))-exact)./max(1,abs(exact)))<2e-6);
end
pass('colored sparse Jacobian versus individual columns');
report.testGroups=numel(report.passed);
if ~exist(fullfile(root,'results'),'dir'),mkdir(fullfile(root,'results'));end
save(fullfile(root,'results','test_report.mat'),'report','-v7');
fprintf('%d test groups passed; scaled RHS error %.3g; runtime %s\n', ...
 report.testGroups,report.rhsError,report.runtime);
 function pass(name),report.passed{end+1}=name;end
end
function reject(f)
failed=false;try,f();catch,failed=true;end
assert(failed,'Expected input rejection.');
end
function [h,S]=exchange(s,~)
h=7*ones(size(s.Tg));S=zeros(size(s.Tg));
end

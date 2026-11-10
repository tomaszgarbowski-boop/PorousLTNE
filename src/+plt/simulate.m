function out = simulate(p,s,o)
%SIMULATE Integrate each constant-control stage using ode15s.
% Time units are seconds. Temperature is Celsius. No optimization toolbox.
% Options: tEnd (86400), target (23), quantile (0.05), stopAtTarget (false),
% maxStep (300 s), initialState ([]), initialStep (0.001 s), useJacobian (true).
p=plt.validate(p);s=plt.schedule(s.times,s.controls);
def=struct('tEnd',86400,'target',23,'quantile',0.05,'stopAtTarget',false, ...
 'maxStep',300,'initialState',[],'initialStep',0.001,'useJacobian',true);
if nargin<3,o=struct();end
assert(isempty(setdiff(fieldnames(o),fieldnames(def))),'plt:options','Unknown option.');
names=fieldnames(def);for k=1:numel(names),if ~isfield(o,names{k}),o.(names{k})=def.(names{k});end,end
assert(o.tEnd>0&&isfinite(o.tEnd)&&o.maxStep>0&&o.initialStep>0&& ...
 o.quantile>=0&&o.quantile<=1&&isfinite(o.target), 'plt:options','Invalid integration options.');
y=plt.initial(p);if ~isempty(o.initialState),y=o.initialState(:);end
assert(numel(y)==5*p.N+7&&all(isfinite(y)),'plt:state','Invalid initial state.');
for k=1:size(s.controls,1)
 [~,rd]=plt.inlet(s.controls(k,:),p);assert(rd>0,'plt:inlet','Nonpositive dry-air density.');
end
checkState(y,p,s.controls(1,:));y0=y;t=[];Y=[];C=[];te=[];ye=[];
if plt.quantile(y(p.N+1:2*p.N),o.quantile)>=o.target,te=0;ye=y';end
edges=[s.times(s.times<o.tEnd);o.tEnd];
timer=tic;
if ~(o.stopAtTarget&&~isempty(te))
 for k=1:numel(edges)-1
  c=s.controls(k,:);
  opts=odeset('RelTol',p.rtol,'AbsTol',p.atol,'JPattern',plt.sparsity(p.N), ...
   'InitialStep',o.initialStep,'MaxStep',o.maxStep,'Events',@(tt,yy) event(tt,yy,p,o,c));
  if o.useJacobian
   opts=odeset(opts,'Jacobian',@(tt,yy) plt.jacobian(tt,yy,c,p));
  end
  [tt,yy,et,ey,ei]=ode15s(@(tt,yy) plt.rhs(tt,yy,c,p),edges(k:k+1),y,opts);
  if any(ei==2),error('plt:domain','Trajectory left the admissible physical domain.');end
  if isempty(te)&&any(ei==1),j=find(ei==1,1);te=et(j);ye=ey(j,:);end
  first=1;if ~isempty(t),first=2;end
  t=[t;tt(first:end)];Y=[Y;yy(first:end,:)];C=[C;repmat(c,numel(tt)-first+1,1)]; %#ok<AGROW>
  y=yy(end,:)';checkState(y,p,c);
  if o.stopAtTarget&&~isempty(te),break;end
 end
else,t=0;Y=y';C=s.controls(1,:);end
out=struct('t',t,'y',Y,'controls',C,'parameters',p,'schedule',s,'options',o, ...
 'eventTime',te,'eventState',ye,'elapsedSeconds',toc(timer),'initialState',y0, ...
 'version','0.2.0','runtime',version);
out.metrics=plt.diagnostics(out);
end
function [v,stop,direction]=event(~,y,p,o,c)
s=plt.stateProperties(y,c,p);
v=[plt.quantile(s.Ts,o.quantile)-o.target;min([s.eg-0.02;0.99-s.eg; ...
 s.X+1e-12;s.qv+1e-12;s.qo+1e-12;s.rd;s.Tg+0.1;60.1-s.Tg;s.Ts+0.1;60.1-s.Ts])];
stop=[double(o.stopAtTarget);1];direction=[1;-1];
end
function checkState(y,p,c)
s=plt.stateProperties(y,c,p);
assert(all(s.eg>0.02&s.eg<0.99)&&all(s.X>=-1e-12)&& ...
 all(s.qv>=-1e-12)&&all(s.qo>=-1e-12)&&all(s.rd>0), ...
 'plt:domain','Invalid gas porosity, concentration, water content or density.');
assert(all(s.Tg>=0&s.Tg<=60)&all(s.Ts>=0&s.Ts<=60), ...
 'plt:domain','The built-in property range is 0 to 60 Celsius.');
end

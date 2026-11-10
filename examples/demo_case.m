function [p,s,o] = demo_case(name,n)
%DEMO_CASE Configurations used by demonstrations and mesh/cost studies.
if nargin<2,n=80;end
switch lower(name)
 case 'heating'
  p=plt.material('inert_granulate');s=plt.schedule(0,[40 .2 .05]);
  o=struct('tEnd',6*3600,'target',35,'maxStep',120);
 case 'phase_cycle'
  p=plt.material('wet_granulate');
  s=plt.schedule([0 4]*3600,[30 .1 .03;40 .98 .03]);
  o=struct('tEnd',8*3600,'target',35,'maxStep',120);
 case 'oxygen_limited'
  p=plt.material('reactive_bed');s=plt.schedule(0,[25 .95 .0001]);
  o=struct('tEnd',12*3600,'target',30,'maxStep',120);
 case 'oxygen_supplied'
  p=plt.material('reactive_bed');s=plt.schedule(0,[25 .95 .0005]);
  o=struct('tEnd',12*3600,'target',30,'maxStep',120);
 otherwise,error('Unknown demonstration.');
end
p.N=n;
end

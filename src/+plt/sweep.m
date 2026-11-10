function runs = sweep(p,s,parameter,values,o)
%SWEEP One-at-a-time study using identical controls and solver options.
% Example: runs=plt.sweep(p,s,'chih',[.2 .4 .8],struct('tEnd',86400));
if nargin<5,o=struct();end
assert(ischar(parameter)&&isfield(p,parameter)&&isnumeric(p.(parameter)), ...
 'plt:sweep','Select an existing numeric scalar parameter.');
assert(isvector(values)&&~isempty(values)&&all(isfinite(values)), ...
 'plt:sweep','Supply a finite nonempty parameter vector.');
runs=cell(numel(values),1);
for i=1:numel(values)
 q=p;q.(parameter)=values(i);runs{i}=plt.simulate(q,s,o);
end
end

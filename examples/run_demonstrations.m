function allRuns = run_demonstrations(n)
%RUN_DEMONSTRATIONS Three physical regimes; reactive case has two flow rates.
if nargin<1,n=80;end
root=fileparts(fileparts(mfilename('fullpath')));
names={'heating','phase_cycle','oxygen_limited','oxygen_supplied'};
allRuns=struct();
for j=1:numel(names)
 name=names{j};[p,s,o]=demo_case(name,n);out=plt.simulate(p,s,o);
 plt.export(out,fullfile(root,'results',name));allRuns.(name)=out;
 fprintf('%s: maxTs=%.4f C, minY=%.5g, wall=%.3f s\n',name, ...
  out.metrics.maxMatrixC,out.metrics.minDryOxygen,out.elapsedSeconds);
end
end

function runs = example_sweep()
%EXAMPLE_SWEEP Sensitivity to interphase heat-transfer correction.
p=plt.defaults();p.N=10;s=plt.schedule(0,[30 .95 .03]);
runs=plt.sweep(p,s,'chih',[.2 .4 .8],struct('tEnd',36*3600));
for i=1:numel(runs),disp(runs{i}.metrics);end
end

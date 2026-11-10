function out = example_staged()
%EXAMPLE_STAGED User-defined policy; illustrative, not an optimized policy.
p=plt.defaults();s=plt.schedule([0 7 14]*3600, ...
 [30 .95 .03;28 .8 .035;26 .65 .04]);
out=plt.simulate(p,s,struct('tEnd',36*3600,'stopAtTarget',true));
disp(out.metrics);plt.export(out,fullfile('results','staged'));
end

function out = example_inert()
%EXAMPLE_INERT Source-free heat transport with inactive phase/reaction terms.
p=plt.defaults();p.bio=false;p.phase=false;
s=plt.schedule(0,[30 0.95 0.03]);
out=plt.simulate(p,s,struct('tEnd',48*3600));disp(out.metrics);
plt.export(out,fullfile('results','inert'));
end

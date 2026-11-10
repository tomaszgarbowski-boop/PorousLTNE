function out = example_compost()
%EXAMPLE_COMPOST Coupled reference bed; Celsius, seconds, SI coefficients.
p=plt.defaults();s=plt.schedule(0,[30 0.95 0.03]);
out=plt.simulate(p,s);disp(out.metrics);
plt.export(out,fullfile('results','compost'));
end

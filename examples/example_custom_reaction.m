function out = example_custom_reaction()
%EXAMPLE_CUSTOM_REACTION Replace biological kinetics without editing solver.
% This deliberately simplified oxygen-limited source is a software example.
p=plt.defaults();p.reactionFcn=@uptake;
out=plt.simulate(p,plt.schedule(0,[30 .95 .03]));disp(out.metrics);
plt.export(out,fullfile('results','custom'));
end
function R = uptake(Ts,X,Y,p)
R=0.5*p.Rbio*exp(-((Ts-25)/12).^2).*X./(0.5+X).*Y./(p.Ko+Y);
end

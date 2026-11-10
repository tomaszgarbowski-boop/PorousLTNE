function out = example_user_exchange()
%EXAMPLE_USER_EXCHANGE A local water-activity and exchange closure.
p=plt.material('wet_granulate');p.N=40;
p.user.aw=.85;p.user.hvol=250;p.user.kmass=.015;
p.exchangeFcn=@user_exchange;
s=plt.schedule([0 2]*3600,[30 .1 .03;40 .98 .03]);
out=plt.simulate(p,s,struct('tEnd',4*3600,'target',35));
root=fileparts(fileparts(mfilename('fullpath')));
plt.export(out,fullfile(root,'results','user_exchange'));
end
function [h,S]=user_exchange(state,p)
h=p.user.hvol*ones(size(state.Ts));
rveq=p.user.aw*plt.psat(state.Ts)./(p.Rv*(state.Ts+273.15));
S=p.user.kmass*(state.rv-rveq);
end

function out = example_user_material()
%EXAMPLE_USER_MATERIAL Custom local kinetics without editing transport code.
p=plt.material('reactive_bed');p.N=40;
p.user.kmax=6e-5; % kg oxygen/(m3 bed s)
p.user.theta=1.06;p.user.Koxygen=.015;p.user.Xhalf=.3;
p.reactionFcn=@user_uptake;
s=plt.schedule(0,[25 .95 .0003]);
out=plt.simulate(p,s,struct('tEnd',6*3600,'target',35));
root=fileparts(fileparts(mfilename('fullpath')));
plt.export(out,fullfile(root,'results','user_material'));
end
function R = user_uptake(T,X,Y,p)
% Local illustrative law: each output entry uses the same cell only.
fT=p.user.theta.^(T-25);fX=X./(p.user.Xhalf+X);
fO=Y./(p.user.Koxygen+Y);
R=p.user.kmax*fT.*fX.*fO;
end

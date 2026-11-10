function y = initial(p)
%INITIAL State blocks Tg, Ts, vapor storage, water/dry mass, oxygen storage.
eg=p.eps-p.rho_bd*p.X0/p.rho_l;
[rv,rd]=plt.inlet([p.T0 p.RH0 0],p);e=ones(p.N,1);
y=[p.T0*e;p.T0*e;eg*rv*e;p.X0*e;eg*rd*p.Yin*e;zeros(7,1)];
end

function [E,W,Es] = inventory(y,p)
%INVENTORY Total sensible/latent enthalpy [MJ], water [kg], wet-solid [MJ].
n=p.N;v=reshape(y(1:5*n),n,5);Tg=v(:,1);Ts=v(:,2);qv=v(:,3);X=v(:,4);
eg=p.eps-p.rho_bd*X/p.rho_l;
rd=(p.p-qv./eg*p.Rv.*(Tg+273.15))./(p.Rda*(Tg+273.15));
Eg=(eg.*rd*p.cda+qv*p.cv).*Tg+qv*p.Lv0;
Es=sum(p.rho_bd*(p.cd+p.cl*X).*Ts)*p.L/n*p.A/1e6;
E=sum(Eg)*p.L/n*p.A/1e6+Es;W=sum(qv+p.rho_bd*X)*p.L/n*p.A;
end

function s = stateProperties(y,c,p)
%PROPERTIES Local constitutive laws and prescribed-flow pressure requirement.
n=p.N;v=reshape(y(1:5*n),n,5);
s.Tg=v(:,1);s.Ts=v(:,2);s.qv=v(:,3);s.X=v(:,4);s.qo=v(:,5);
s.eg=p.eps-p.rho_bd*s.X/p.rho_l;eg0=p.eps-p.rho_bd*p.X0/p.rho_l;
s.rv=s.qv./s.eg;
s.rd=(p.p-s.rv.*p.Rv.*(s.Tg+273.15))./(p.Rda*(s.Tg+273.15));
s.rg=s.rd+s.rv;[s.vin,s.din,s.gda]=plt.inlet(c,p);
s.u=s.gda./s.rd;s.Y=s.qo./(s.eg.*s.rd);
s.mu=p.mu0*((s.Tg+273.15)/293.15).^0.7;
s.K=p.K0*(s.eg/eg0).^3.*((1-eg0)./(1-s.eg)).^2;
Re=s.rg.*s.u*p.deq./(s.mu.*s.eg);Nu=2+1.1*Re.^0.6*p.Pr^(1/3);
massSpeed=(2*p.Dv+1.1*Re.^0.6.*(s.mu./s.rg).^(1/3)*p.Dv^(2/3))/p.deq;
area=p.chia*6*(1-p.eps)/p.deq;
s.hvol=p.chih*Nu*p.kair/p.deq*area;
s.Sp=p.chim*massSpeed*area.* ...
 (s.rv-plt.psat(s.Ts)./(p.Rv*(s.Ts+273.15)));
if ~p.phase,s.Sp=zeros(n,1);end
if ~isempty(p.exchangeFcn)
 [s.hvol,s.Sp]=p.exchangeFcn(s,p);
 assert(isequal(size(s.hvol),[n 1])&&isequal(size(s.Sp),[n 1])&& ...
  all(isfinite(s.hvol))&&all(s.hvol>=0)&&all(isfinite(s.Sp)), ...
  'plt:exchange','Exchange callback must return two finite N-by-1 columns.');
 if ~p.phase,s.Sp=zeros(n,1);end
end
Y=s.Y;if ~p.oxygen,Y=p.Yin*ones(n,1);end
width=p.sigp*ones(n,1);width(s.Ts<p.Topt)=p.sigm;
s.R=p.Rbio*exp(-((s.Ts-p.Topt)./width).^2).* ...
 exp(-((s.X-p.Xopt)/p.sigX).^2).*Y./(p.Ko+Y);
if ~isempty(p.reactionFcn),s.R=p.reactionFcn(s.Ts,s.X,Y,p);end
if ~p.bio,s.R=zeros(n,1);end
assert(isequal(size(s.R),[n 1])&&all(isfinite(s.R))&&all(s.R>=0), ...
 'plt:reaction','Reaction callback must return finite nonnegative N-by-1 uptake.');
s.dp=sum(s.mu./s.K.*s.u+s.rg*p.CF./sqrt(s.K).*s.u.^2)*(p.L/n);
end

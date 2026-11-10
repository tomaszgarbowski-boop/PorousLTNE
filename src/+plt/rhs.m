function [dy,localRates] = rhs(~,y,c,p)
%RHS Five field blocks plus seven integrated boundary/source diagnostics.
% Conservative gas/wet-solid enthalpy, including transported water enthalpy.
s=plt.stateProperties(y,c,p);n=p.N;dz=p.L/n;
Tf=[c(1);s.Tg];hv=p.Lv0+p.cv*Tf;
Fva=s.gda*[s.vin/s.din;s.rv./s.rd];
Jv=plt.diffFlux(s.rv,p.Dv*s.eg.^(1+p.nD),dz,s.vin);Fv=Fva+Jv;
dqv=-diff(Fv)/dz-s.Sp;
Jw=plt.diffFlux(s.X,p.rho_bd*p.Dw,dz);
dX=(-diff(Jw)/dz+s.Sp)/p.rho_bd;
Fo=s.gda*[p.Yin;s.Y]+plt.diffFlux(s.Y,p.Do*s.rd.*s.eg.^(1+p.nO),dz,p.Yin);
dqo=-diff(Fo)/dz-s.R;if ~p.oxygen,dqo=zeros(n,1);end
Td=[c(1);0.5*(s.Tg(1:end-1)+s.Tg(2:end));s.Tg(end)];
Fg=s.gda*p.cda*Tf+Fva.*hv+Jv.*(p.Lv0+p.cv*Td)+plt.diffFlux(s.Tg,p.kg,dz,c(1));
Hpc=s.Sp.*(p.Lv0+p.cv*s.Ts);Qgs=s.hvol.*(s.Tg-s.Ts);
rg=-diff(Fg)/dz-Qgs-Hpc;deg=-p.rho_bd/p.rho_l*dX;
B=p.cv-p.Rv/p.Rda*p.cda;
ET=s.eg*p.p*p.cda*273.15./(p.Rda*(s.Tg+273.15).^2)+s.qv*B;
Eq=B*s.Tg+p.Lv0;Ee=p.p*p.cda/p.Rda*s.Tg./(s.Tg+273.15);
dTg=(rg-Eq.*dqv-Ee.*deg)./ET;
Tw=[s.Ts(1);0.5*(s.Ts(1:end-1)+s.Ts(2:end));s.Ts(end)];
Fs=plt.diffFlux(s.Ts,p.ks*(1+p.alphak*(s.X-p.X0)),dz)+Jw.*p.cl.*Tw;
Qwall=p.Uw*p.AwV*(p.Tamb-s.Ts);Qbio=p.Ho*s.R;
rs=-diff(Fs)/dz+Qgs+Hpc+Qbio+Qwall;
Cs=p.rho_bd*(p.cd+p.cl*s.X);
dTs=(rs-p.rho_bd*p.cl*s.Ts.*dX)./Cs;
[va,da]=plt.inlet([p.Tamb p.RHamb c(3)],p);
hin=p.cda*c(1)+s.vin/s.din*(p.Lv0+p.cv*c(1));
hamb=p.cda*p.Tamb+va/da*(p.Lv0+p.cv*p.Tamb);
% Each diagnostic contribution is assigned to the cell on which it depends.
% This also permits colored differentiation of the global diagnostic rows.
localRates=zeros(7,n);
localRates(1,1)=Fg(1)*p.A/1e6;
localRates(1,n)=localRates(1,n)-Fg(end)*p.A/1e6;
localRates(2,:)=Qbio'*dz*p.A/1e6;
localRates(3,:)=Qwall'*dz*p.A/1e6;
localRates(4,1)=Fv(1)*p.A;
localRates(4,n)=localRates(4,n)-Fv(end)*p.A;
localRates(5,:)=max(s.Sp,0)'*dz*p.A;
localRates(6,:)=(max(s.Sp,0).*(p.Lv0+(p.cv-p.cl)*s.Ts))'*dz*p.A/1e6;
dpi=(s.mu./s.K.*s.u+s.rg*p.CF./sqrt(s.K).*s.u.^2)*dz;
localRates(7,:)=dpi'*c(3)*p.A/p.eta/3.6e6;
localRates(7,1)=localRates(7,1)+max(0,s.gda*p.A*(hin-hamb))/3.6e6;
rates=sum(localRates,2);
dy=[dTg;dTs;dqv;dX;dqo;rates];
end

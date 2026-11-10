"""Conservative 1D LTNE compost model, audited 2026-11-10.
Enthalpy reference: dry air/liquid/dry solid at 0 degC; vapor h_v(0)=2.501 MJ/kg.
Oxygen is transported per kg dry carrier gas. No experimental calibration.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

@dataclass
class Params:
    N:int=40
    L:float=1.; A:float=1.; eps:float=.8875; rho_bd:float=180.; rho_l:float=997.
    p:float=101325.; Rda:float=287.05; Rv:float=461.5
    cda:float=1006.; cv:float=1860.; cl:float=4180.; cd:float=1500.; Lv0:float=2.501e6
    X0:float=2.; T0:float=3.; RH0:float=.95; Tamb:float=23.; RHamb:float=.5
    K0:float=2e-9; CF:float=.55; deq:float=.03; mu0:float=1.8e-5
    ks:float=.2; alphak:float=.25; kg:float=.03; kair:float=.026; Pr:float=.71
    chia:float=.3; chih:float=.4; chim:float=.02; Dv:float=2.4e-5; Do:float=2e-5
    nD:float=1.5; nO:float=1.5; Dw:float=2e-10
    Rbio:float=2e-5; Topt:float=25.; sigm:float=10.; sigp:float=8.; Xopt:float=2.1; sigX:float=.8
    Ko:float=.02; Ho:float=13.6e6; Yin:float=.233; Uw:float=.3; AwV:float=4.; eta:float=.6
    phase:bool=True; bio:bool=True; oxygen:bool=True
    rtol:float=2e-6; atol:float=1e-9
    method:str='BDF'

def psat(T):
    return 611.21*np.exp((18.678-T/234.5)*T/(257.14+T))

class Model:
    def __init__(self,p=Params()):
        self.p=p;self.N=p.N;self.dz=p.L/p.N
        self.eps0=p.eps-p.rho_bd*p.X0/p.rho_l
        n=5*p.N+7; sp=lil_matrix((n,n),dtype=int)
        for i in range(p.N):
            for a in range(5):
                for b in range(5):
                    for j in range(max(0,i-1),min(p.N,i+2)):sp[a*p.N+i,b*p.N+j]=1
        sp[5*p.N:,:5*p.N]=1
        self.sparsity=sp.tocsr()
    def inlet(self,c):
        p=self.p;T,RH,u=c;v=RH*psat(T)/(p.Rv*(T+273.15));da=(p.p-v*p.Rv*(T+273.15))/(p.Rda*(T+273.15))
        return v,da,da*u
    def initial(self):
        p=self.p;v,da,_=self.inlet((p.T0,p.RH0,.03))
        return np.r_[np.full(self.N,p.T0),np.full(self.N,p.T0),np.full(self.N,self.eps0*v),np.full(self.N,p.X0),np.full(self.N,self.eps0*da*p.Yin),np.zeros(7)]
    def props(self,y,c):
        p=self.p;N=self.N;Tg,Ts,qv,X,qo=y[:5*N].reshape(5,N);eg=p.eps-p.rho_bd*X/p.rho_l
        rv=qv/eg;rd=(p.p-rv*p.Rv*(Tg+273.15))/(p.Rda*(Tg+273.15));rg=rd+rv
        vin,din,gda=self.inlet(c);u=gda/rd;Y=qo/(eg*rd)
        mu=p.mu0*((Tg+273.15)/293.15)**.7
        K=p.K0*(eg/self.eps0)**3*((1-self.eps0)/(1-eg))**2
        Re=rg*u*p.deq/(mu*eg);Nu=2+1.1*Re**.6*p.Pr**(1/3);Sc=mu/(rg*p.Dv);Sh=2+1.1*Re**.6*Sc**(1/3)
        area=p.chia*6*(1-p.eps)/p.deq;hvol=p.chih*Nu*p.kair/p.deq*area
        Sp=p.chim*Sh*p.Dv/p.deq*area*(rv-psat(Ts)/(p.Rv*(Ts+273.15))) if p.phase else np.zeros(N)
        fT=np.exp(-((Ts-p.Topt)/np.where(Ts<p.Topt,p.sigm,p.sigp))**2)
        Yb=Y if p.oxygen else p.Yin
        R=p.Rbio*fT*np.exp(-((X-p.Xopt)/p.sigX)**2)*Yb/(p.Ko+Yb) if p.bio else np.zeros(N)
        dp=(mu/K*u+rg*p.CF/np.sqrt(K)*u*u).sum()*self.dz
        return Tg,Ts,qv,X,qo,eg,rv,rd,rg,Y,u,vin,din,gda,hvol,Sp,R,K,dp
    def diff_flux(self,x,k,left=None):
        k=np.broadcast_to(k,(self.N,));f=np.zeros(self.N+1)
        f[1:-1]=-2*k[:-1]*k[1:]/(k[:-1]+k[1:])*(x[1:]-x[:-1])/self.dz
        if left is not None:f[0]=-k[0]*(x[0]-left)/(.5*self.dz)
        return f
    def rhs(self,t,y,c):
        p=self.p;N=self.N
        Tg,Ts,qv,X,qo,eg,rv,rd,rg,Y,u,vin,din,gda,hvol,Sp,R,K,dp=self.props(y,c)
        Tface=np.r_[c[0],Tg];hv=p.Lv0+p.cv*Tface
        Fv_adv=gda*np.r_[vin/din,rv/rd]
        Jv=self.diff_flux(rv,eg*p.Dv*eg**p.nD,vin)
        Fv=Fv_adv+Jv;dqv=-np.diff(Fv)/self.dz-Sp
        Jw=self.diff_flux(X,p.rho_bd*p.Dw)
        dX=(-np.diff(Jw)/self.dz+Sp)/p.rho_bd
        Fo=gda*np.r_[p.Yin,Y]+self.diff_flux(Y,eg*rd*p.Do*eg**p.nO,p.Yin)
        dqo=-np.diff(Fo)/self.dz-R if p.oxygen else np.zeros(N)
        # Total gas enthalpy flux: advection, vapor diffusion, heat conduction.
        Tdiff=np.r_[c[0],.5*(Tg[:-1]+Tg[1:]),Tg[-1]]
        Fg=gda*p.cda*Tface+Fv_adv*hv+Jv*(p.Lv0+p.cv*Tdiff)+self.diff_flux(Tg,p.kg,c[0])
        # Interphase mass carries vapor enthalpy at the wet-matrix temperature.
        Hpc=Sp*(p.Lv0+p.cv*Ts);Qgs=hvol*(Tg-Ts)
        rhs_g=-np.diff(Fg)/self.dz-Qgs-Hpc
        deg=-p.rho_bd/p.rho_l*dX
        B=p.cv-p.Rv/p.Rda*p.cda
        dEdT=eg*p.p*p.cda*273.15/(p.Rda*(Tg+273.15)**2)+qv*B
        dEdq=B*Tg+p.Lv0;dEde=p.p*p.cda/p.Rda*Tg/(Tg+273.15)
        dTg=(rhs_g-dEdq*dqv-dEde*deg)/dEdT
        Twface=np.r_[Ts[0],.5*(Ts[:-1]+Ts[1:]),Ts[-1]]
        Fs=self.diff_flux(Ts,p.ks*(1+p.alphak*(X-p.X0)))+Jw*p.cl*Twface
        Qwall=p.Uw*p.AwV*(p.Tamb-Ts);Qbio=p.Ho*R
        rhs_s=-np.diff(Fs)/self.dz+Qgs+Hpc+Qbio+Qwall
        Cs=p.rho_bd*(p.cd+p.cl*X)
        dTs=(rhs_s-p.rho_bd*p.cl*Ts*dX)/Cs
        va,da,_=self.inlet((p.Tamb,p.RHamb,c[2]));hin=p.cda*c[0]+vin/din*(p.Lv0+p.cv*c[0]);hamb=p.cda*p.Tamb+va/da*(p.Lv0+p.cv*p.Tamb)
        # Integrals: net boundary enthalpy MJ, bio MJ, wall MJ, net water kg,
        # positive condensation kg, positive latent heat MJ, conditioning kWh.
        rates=[(Fg[0]-Fg[-1])*p.A/1e6,Qbio.sum()*self.dz*p.A/1e6,Qwall.sum()*self.dz*p.A/1e6,(Fv[0]-Fv[-1])*p.A,np.maximum(Sp,0).sum()*self.dz*p.A,(np.maximum(Sp,0)*(p.Lv0+(p.cv-p.cl)*Ts)).sum()*self.dz*p.A/1e6,(max(0,gda*p.A*(hin-hamb))+dp*c[2]*p.A/p.eta)/3.6e6]
        return np.r_[dTg,dTs,dqv,dX,dqo,rates]
    def energy(self,y):
        p=self.p;N=self.N;Tg,Ts,qv,X,qo=y[:5*N].reshape(5,N);eg=p.eps-p.rho_bd*X/p.rho_l
        rd=(p.p-qv/eg*p.Rv*(Tg+273.15))/(p.Rda*(Tg+273.15))
        Eg=(eg*rd*p.cda+qv*p.cv)*Tg+qv*p.Lv0
        Es=p.rho_bd*(p.cd+p.cl*X)*Ts
        return (Eg.sum()+Es.sum())*self.dz*p.A/1e6,Es.sum()*self.dz*p.A/1e6,(qv.sum()+p.rho_bd*X.sum())*self.dz*p.A
    def run(self,controls=None,horizon=48.,stop=False,sample=.1):
        if controls is None:controls=[(30.,.95,.03)]*3
        if len(controls)==1:controls=controls*3
        p=self.p;N=self.N;y0=self.initial();y=y0.copy();ts=[];ys=[];cs=[];t95=None;yevent=None
        bounds=[0.,7.,14.,horizon];bounds=sorted(set(x for x in bounds if x<=horizon))
        for a,b in zip(bounds[:-1],bounds[1:]):
            c=controls[min(int(a>=7)+int(a>=14),2)]
            def event(t,y):return np.quantile(y[N:2*N],.05,method='linear')-23
            event.direction=1;event.terminal=stop
            sol=solve_ivp(lambda t,y:self.rhs(t,y,c),(a*3600,b*3600),y,method=p.method,rtol=p.rtol,atol=p.atol,jac_sparsity=self.sparsity,events=event,dense_output=True)
            if not sol.success:raise RuntimeError(sol.message)
            if len(sol.t_events[0]) and t95 is None:t95=sol.t_events[0][0]/3600;yevent=sol.y_events[0][0]
            tt=np.unique(np.r_[np.arange(a*3600,sol.t[-1],sample*3600),sol.t,sol.t_events[0]])
            if ts and tt[0]==ts[-1]:tt=tt[1:]
            yy=sol.sol(tt);ts.extend(tt/3600);ys.extend(yy.T);cs.extend([c]*len(tt));y=sol.y[:,-1]
            if stop and t95 is not None:break
        ts=np.array(ts);ys=np.array(ys);cs=np.array(cs)
        if yevent is None:yevent=ys[-1]
        idx=ts<=(t95 if t95 is not None else ts[-1])+1e-8
        metrics={'t95':t95,'E':float(yevent[-1]),'Tmax':float(ys[idx,N:2*N].max()),'Xmin':float(ys[idx,3*N:4*N].min()),'cond':float(yevent[-3]),'latent':float(yevent[-2]),'bio':float(yevent[-6]),'stored':self.energy(yevent)[1]-self.energy(y0)[1]}
        props=[self.props(v,c) for v,c in zip(ys,cs)]
        metrics['dpmax']=max(s[-1] for s,yes in zip(props,idx) if yes)
        metrics['Ymin']=min(s[9].min() for s,yes in zip(props,idx) if yes)
        metrics['Ywetmin']=min((s[9]*s[7]/s[8]).min() for s,yes in zip(props,idx) if yes)
        e0=self.energy(y0);ef=self.energy(ys[-1]);ints=ys[-1,-7:]
        metrics['water_err']=100*abs(ef[2]-e0[2]-ints[3])/max(abs(ints[3]),1e-8)
        metrics['energy_err']=100*abs(ef[0]-e0[0]-sum(ints[:3]))/max(abs(ints[0])+abs(ints[1])+abs(ints[2]),1e-8)
        metrics['dTmax']=float(np.max(abs(ys[idx,:N]-ys[idx,N:2*N])))
        metrics['dTmean']=float(np.trapezoid(np.mean(abs(ys[idx,:N]-ys[idx,N:2*N]),axis=1),ts[idx])/ts[idx][-1])
        return {'t':ts,'y':ys,'c':cs,'metrics':metrics,'event':yevent,'N':N}

if __name__=='__main__':
    import time,json
    for n in [10,40]:
        s=time.time();r=Model(Params(N=n)).run(horizon=24);print(n,time.time()-s,json.dumps(r['metrics']),flush=True)

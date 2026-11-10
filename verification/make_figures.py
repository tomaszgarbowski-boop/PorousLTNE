from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
ROOT=Path(__file__).resolve().parents[1]
PKG=ROOT;DEST=PKG/'figures';DEST.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
colors=['#2466A4','#D17B23','#2E8B72']
def save(fig,name):
 for ext in ['png','svg','pdf']:fig.savefig(DEST/f'{name}.{ext}',dpi=400,bbox_inches='tight',facecolor='white')
 plt.close(fig)
def panels(axs):
 for i,a in enumerate(np.asarray(axs).ravel()):
  a.text(-.15,1.04,f'({chr(97+i)})',transform=a.transAxes,fontweight='bold');a.grid(alpha=.2)
def get(name):return loadmat(PKG/f'results/{name}/result.mat',simplify_cells=True)['out']
def fields(o):
 p=o['parameters'];N=int(p['N']);y=o['y'];Tg=y[:,:N];Ts=y[:,N:2*N];qv=y[:,2*N:3*N];X=y[:,3*N:4*N];qo=y[:,4*N:5*N]
 eg=p['eps']-p['rho_bd']*X/p['rho_l'];rv=qv/eg;rd=(p['p']-rv*p['Rv']*(Tg+273.15))/(p['Rda']*(Tg+273.15));Y=qo/(eg*rd)
 W=(qv+p['rho_bd']*X).sum(axis=1)*p['A']*p['L']/N
 return p,N,o['t']/3600,Tg,Ts,X,Y,W
fig,ax=plt.subplots(figsize=(6.5,3.2));ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off')
labels=[(.2,4,'Material and state\nproperties and closures\nlocal user constants'),(3.65,4,'Inlet schedule\ntemperature and RH\nprescribed airflow'),(7.2,4,'Conservative core\nfive coupled fields\nenthalpy and species'),(7.2,1.5,'Stage integration\nsparse Jacobian\node15s and events'),(3.65,1.5,'Diagnostics\nfields and balances\nMAT and CSV export'),(.2,1.5,'Verification and reuse\nexact benchmarks\nthree demonstrations')]
for x,y,txt in labels:
 ax.add_patch(FancyBboxPatch((x,y),2.5,1.3,boxstyle='round,pad=.08',fc='#EDF3F8',ec='#38566F'))
 ax.text(x+1.25,y+.65,txt,ha='center',va='center',fontsize=7.3)
for a,b in [((2.8,4.65),(3.55,4.65)),((6.25,4.65),(7.1,4.65)),((8.45,3.9),(8.45,2.9)),((7.1,2.15),(6.3,2.15)),((3.55,2.15),(2.8,2.15))]:
 ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,color='#38566F',lw=1.2))
ax.text(5,.35,'Shared transport equations; replaceable local material, reaction and exchange laws',ha='center',fontsize=7.5)
save(fig,'fig1_architecture')
b=loadmat(PKG/'results/benchmark_report.mat',simplify_cells=True)['report']
fig,axs=plt.subplots(2,2,figsize=(6.5,4.8),layout='constrained')
d=b['diffusion'];a=b['advection'];r=b['reaction']
axs[0,0].loglog(d[:,0],d[:,1],'o-',color=colors[0],label='RMS error')
axs[0,0].loglog(d[:,0],d[0,1]*(d[0,0]/d[:,0])**2,'--',color='.5',label=r'$N^{-2}$')
axs[0,0].set(xlabel='Cells, N',ylabel='Temperature error [°C]');axs[0,0].legend(frameon=False)
axs[0,1].loglog(a[:,0],a[:,1],'o-',color=colors[1],label='L1 error')
axs[0,1].loglog(a[:,0],a[-1,1]*(a[-1,0]/a[:,0]),'--',color='.5',label=r'$N^{-1}$')
axs[0,1].set(xlabel='Cells, N',ylabel='Oxygen-ratio error [kg/kg]');axs[0,1].legend(frameon=False)
axs[1,0].loglog(r[:,0],r[:,1],'o-',color=colors[2]);axs[1,0].set(xlabel='Relative tolerance',ylabel='Normalized decay error')
f=b['advectionProfile'];axs[1,1].plot(f[:,0],f[:,1],color='.2',label='Exact cell average');axs[1,1].plot(f[:,0],f[:,2],'--',color=colors[1],label='Numerical, N = 160')
axs[1,1].set(xlabel='Axial position [m]',ylabel='Oxygen ratio [kg/kg]');axs[1,1].legend(frameon=False,loc='upper left',fontsize=6.2)
panels(axs);save(fig,'fig2_benchmarks')
fig,axs=plt.subplots(1,2,figsize=(6.5,2.8),layout='constrained');o=get('heating');p,N,t,Tg,Ts,X,Y,W=fields(o)
axs[0].plot(t,Tg.mean(axis=1),color=colors[1],label='Mean gas')
axs[0].plot(t,Ts.mean(axis=1),color=colors[0],label='Mean matrix')
axs[0].fill_between(t,Ts.min(axis=1),Ts.max(axis=1),color=colors[0],alpha=.15,label='Matrix range')
axs[0].set(xlabel='Time [h]',ylabel='Temperature [°C]',xlim=(0,6));axs[0].legend(frameon=False)
z=(np.arange(N)+.5)*p['L']/N
for tt,col in zip([.25,.5,1,2,6],plt.cm.viridis(np.linspace(.1,.9,5))):
 yy=np.array([np.interp(tt,t,Ts[:,i]) for i in range(N)]);axs[1].plot(z,yy,label=f'{tt:g} h',color=col)
axs[1].set(xlabel='Axial position [m]',ylabel='Matrix temperature [°C]');axs[1].legend(frameon=False,fontsize=6.5)
panels(axs);save(fig,'fig3_heating')
fig,axs=plt.subplots(1,2,figsize=(6.5,2.8),layout='constrained');o=get('phase_cycle');p,N,t,Tg,Ts,X,Y,W=fields(o)
axs[0].plot(t,Tg.mean(axis=1),color=colors[1],label='Mean gas');axs[0].plot(t,Ts.mean(axis=1),color=colors[0],label='Mean matrix');axs[0].step([0,4,8],[30,40,40],where='post',ls='--',color='.5',label='Inlet')
axs[0].set(xlabel='Time [h]',ylabel='Temperature [°C]');axs[0].legend(frameon=False)
axs[1].plot(t,W-W[0],color=colors[2],label='Total stored water')
wl=X.sum(axis=1)*p['rho_bd']*p['A']*p['L']/N
axs[1].plot(t,wl-wl[0],'--',color=colors[1],label='Liquid water')
axs[1].set(xlabel='Time [h]',ylabel='Water-inventory change [kg]');axs[1].legend(frameon=False)
for a in axs:a.axvline(4,color='.5',ls=':',lw=1);a.set_xlim(0,8)
panels(axs);save(fig,'fig4_phase_cycle')
fig,axs=plt.subplots(1,2,figsize=(6.5,2.8),layout='constrained')
for name,col,label in [('oxygen_limited',colors[0],'u = 0.1 mm/s'),('oxygen_supplied',colors[1],'u = 0.5 mm/s')]:
 p,N,t,Tg,Ts,X,Y,W=fields(get(name));axs[0].plot(t,Ts.mean(axis=1),color=col,label=label)
 axs[1].plot(t,Y[:,-1]/p['Yin'],color=col,label=label)
axs[0].set(xlabel='Time [h]',ylabel='Mean matrix temperature [°C]');axs[1].set(xlabel='Time [h]',ylabel='Outlet / inlet oxygen ratio')
for a in axs:a.set_xlim(0,12);a.legend(frameon=False)
panels(axs);save(fig,'fig5_oxygen')
m=loadmat(PKG/'results/mesh_cost_report.mat',simplify_cells=True)['study']
fig,axs=plt.subplots(1,2,figsize=(6.5,2.8),layout='constrained')
for name,col,label in zip(['heating','phase_cycle','oxygen_limited'],colors,['Heating','Phase cycle','Reactive bed']):
 a=m[name]['rows'];axs[0].errorbar(a[:,0],a[:,4],yerr=[a[:,4]-a[:,5],a[:,6]-a[:,4]],marker='o',color=col,label=label,capsize=3)
 axs[1].loglog(a[:,0],a[:,1],'o-',color=col,label=label)
axs[0].set(xlabel='Cells, N',ylabel='Median wall time [s]');axs[0].set_xscale('log',base=2);axs[0].set_yscale('log')
axs[1].set(xlabel='Cells, N',ylabel='Maximum matrix-T difference [°C]')
for a in axs:a.legend(frameon=False);a.set_xticks([10,20,40,80]);a.set_xticklabels(['10','20','40','80']);a.xaxis.set_minor_formatter(NullFormatter())
panels(axs);save(fig,'fig6_cost_accuracy')
print('Created six figure sets in',DEST)

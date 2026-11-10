function p = material(name)
%MATERIAL Reproducible illustrative presets, not measured material data.
% Names: compost, inert_granulate, wet_granulate, reactive_bed.
p=plt.defaults();
switch lower(name)
 case 'compost'
 case 'inert_granulate'
  p.L=.4;p.rho_bd=450;p.eps=.65;p.X0=.15;p.T0=15;
  p.cd=1100;p.ks=.35;p.alphak=0;p.deq=.008;p.K0=8e-9;
  p.chia=.6;p.chih=.5;p.bio=false;p.phase=false;p.Uw=0;
 case 'wet_granulate'
  p.L=.3;p.rho_bd=250;p.eps=.7;p.X0=.8;p.T0=25;p.RH0=1;
  p.cd=1200;p.ks=.25;p.alphak=0;p.deq=.01;p.K0=6e-9;
  p.chia=.6;p.chih=.6;p.chim=.05;p.bio=false;p.Uw=0;
 case 'reactive_bed'
  p.L=.6;p.rho_bd=120;p.eps=.85;p.X0=1;p.T0=25;
  p.Tamb=25;p.cd=1400;p.ks=.12;p.alphak=0;p.deq=.015;
  p.Topt=35;p.sigm=15;p.sigp=8;p.Xopt=1;p.sigX=.6;
  p.Rbio=8e-5;p.Yin=.08;p.Ko=.02;p.Uw=.15;
 otherwise,error('plt:material','Unknown preset: %s',name);
end
p=plt.validate(p);
end

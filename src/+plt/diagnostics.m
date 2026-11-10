function m = diagnostics(out)
%DIAGNOSTICS Event metrics and global balances over the simulated horizon.
p=out.parameters;n=p.N;Y=out.y;t=out.t;
if isempty(out.eventTime),m.activationHours=NaN;v=Y(end,:)';else
 m.activationHours=out.eventTime/3600;v=out.eventState';end
m.conditioningKWh=v(end);m.condensedKg=v(end-2);m.latentMJ=v(end-1);
m.bioMJ=v(end-5);
[E0,W0]=plt.inventory(out.initialState,p);[Ef,Wf]=plt.inventory(Y(end,:)',p);
I=Y(end,end-6:end);
m.energyResidualMJ=Ef-E0-sum(I(1:3));m.waterResidualKg=Wf-W0-I(4);
m.energyResidualPercent=100*abs(m.energyResidualMJ)/max(sum(abs(I(1:3))),1e-8);
m.waterResidualPercent=100*abs(m.waterResidualKg)/max(abs(I(4)),1e-8);
m.maxMatrixC=max(max(Y(:,n+1:2*n)));m.minWater=min(min(Y(:,3*n+1:4*n)));
m.maxPressurePa=0;m.minDryOxygen=Inf;m.maxPhaseGapC=max(max(abs(Y(:,1:n)-Y(:,n+1:2*n))));
for i=1:numel(t)
 s=plt.stateProperties(Y(i,:)',out.controls(i,:),p);
 m.maxPressurePa=max(m.maxPressurePa,s.dp);m.minDryOxygen=min(m.minDryOxygen,min(s.Y));
end
end

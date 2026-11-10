function h = plotResult(out)
%PLOTRESULT Standard temperature and water/energy diagnostic plots.
% Requires a graphics-capable MATLAB or Octave session.
n=out.parameters.N;t=out.t/3600;Ts=out.y(:,n+1:2*n);
h=figure('Color','w');
subplot(2,2,1);plot(t,[min(Ts,[],2),mean(Ts,2),max(Ts,[],2)],'LineWidth',1.3);
xlabel('Time [h]');ylabel('Matrix temperature [C]');legend('Minimum','Mean','Maximum','Location','best');grid on;
subplot(2,2,2);plot(((1:n)-.5)/n,Ts(end,:),'LineWidth',1.3);
xlabel('Normalized position z/L');ylabel('Final matrix temperature [C]');grid on;
subplot(2,2,3);plot(t,out.y(:,end-2),'LineWidth',1.3);
xlabel('Time [h]');ylabel('Cumulative condensation [kg]');grid on;
subplot(2,2,4);plot(t,out.y(:,end),'LineWidth',1.3);
xlabel('Time [h]');ylabel('Conditioning indicator [kWh]');grid on;
end

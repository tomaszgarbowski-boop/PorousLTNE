function [rv,rd,gda] = inlet(c,p)
%INLET c = [temperature_C, relative_humidity_fraction, velocity_m_per_s].
rv=c(2)*plt.psat(c(1))/(p.Rv*(c(1)+273.15));
rd=(p.p-rv*p.Rv*(c(1)+273.15))/(p.Rda*(c(1)+273.15));
gda=rd*c(3);
end

function f = diffFlux(x,k,dz,left)
%DIFFFLUX Conservative face flux with harmonic internal coefficients.
% Closed right boundary; optional Dirichlet left boundary at half-cell distance.
n=numel(x);if isscalar(k),k=repmat(k,n,1);end
f=zeros(n+1,1);den=k(1:end-1)+k(2:end);kh=zeros(n-1,1);
active=den>0;kh(active)=2*k(find(active)).*k(find(active)+1)./den(active);
f(2:n)=-kh.*diff(x)/dz;
if nargin==4,f(1)=-k(1)*(x(1)-left)/(0.5*dz);end
end

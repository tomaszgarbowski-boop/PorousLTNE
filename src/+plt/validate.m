function p = validate(p)
%VALIDATE Reject incomplete, misspelled and physically invalid parameters.
q = plt.defaults(); expected = fieldnames(q); actual = fieldnames(p);
assert(isempty(setxor(expected,actual)), 'plt:parameters', ...
 'Use plt.defaults() and modify known fields only.');
assert(isscalar(p.N) && isfinite(p.N) && p.N>=2 && p.N==fix(p.N), ...
 'plt:grid','N must be an integer >= 2.');
positive = {'L','A','rho_bd','rho_l','p','Rda','Rv','cda','cv','cl','cd', ...
 'Lv0','K0','deq','mu0','kair','Pr','sigm','sigp', ...
 'sigX','Ko','eta','rtol','atol'};
for i=1:numel(positive)
 v=p.(positive{i}); assert(isscalar(v)&&isfinite(v)&&v>0,'plt:positive', ...
  '%s must be a finite positive scalar.',positive{i});
end
nonnegative={'X0','CF','alphak','chia','chih','chim','nD','nO','Dw', ...
 'Rbio','Xopt','Yin','Uw','AwV','ks','kg','Dv','Do','Ho'};
for i=1:numel(nonnegative)
 v=p.(nonnegative{i}); assert(isscalar(v)&&isfinite(v)&&v>=0, ...
  'plt:nonnegative','%s must be a finite nonnegative scalar.',nonnegative{i});
end
for name={'T0','Tamb','Topt'}
 v=p.(name{1});assert(isscalar(v)&&isfinite(v)&&v>=0&&v<=60, ...
  'plt:temperature','Reference temperature must be in [0,60] Celsius.');
end
assert(p.eps>0&&p.eps<1&&p.eps-p.rho_bd*p.X0/p.rho_l>0.02, ...
 'plt:porosity','Initial gas porosity must exceed 0.02 and total porosity be <1.');
assert(p.RH0>=0&&p.RH0<=1&&p.RHamb>=0&&p.RHamb<=1&&p.eta<=1, ...
 'plt:fraction','Relative humidity and fan efficiency must not exceed one.');
for name={'phase','bio','oxygen'}
 v=p.(name{1});assert(isscalar(v)&&(v==0||v==1),'plt:switch','Invalid process switch.');
end
for name={'reactionFcn','exchangeFcn'}
 assert(isempty(p.(name{1}))||isa(p.(name{1}),'function_handle'), ...
  'plt:callback','Callback must be empty or a function handle.');
end
assert(isstruct(p.user)&&isscalar(p.user),'plt:user','user must be a scalar struct.');
for name=fieldnames(p.user)'
 v=p.user.(name{1});
 assert((isnumeric(v)||islogical(v))&&isreal(v)&&all(isfinite(v(:))), ...
  'plt:user','User constants must be finite real numeric arrays.');
end
[~,rd]=plt.inlet([p.T0 p.RH0 0],p);assert(rd>0,'plt:pressure','Invalid gas state.');
end

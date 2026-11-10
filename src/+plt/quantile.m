function q = quantile(x,a)
%QUANTILE Explicit linear order-statistic interpolation (type 7).
assert(a>=0&&a<=1&&~isempty(x),'plt:quantile','Invalid quantile.');
x=sort(x(:));h=1+(numel(x)-1)*a;k=floor(h);j=min(k+1,numel(x));
q=x(k)+(h-k)*(x(j)-x(k));
end

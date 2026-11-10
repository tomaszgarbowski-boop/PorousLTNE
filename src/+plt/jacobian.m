function J = jacobian(t,y,c,p)
%JACOBIAN Sparse forward-difference Jacobian using 15 structural colors.
% Five local fields and three spatial colors cover nearest-neighbour coupling.
% Seven diagnostic rows are recovered from cell-local integral contributions.
n=p.N;m=5*n+7;[f,g]=plt.rhs(t,y,c,p);
scale=kron([1;1;1e-3;1;1e-2],ones(n,1));
h=sqrt(eps)*max(abs(y(1:5*n)),scale);
ri=zeros(110*n,1);ci=ri;vv=ri;cursor=0;
for field=0:4
 for color=1:3
  cells=color:3:n;cols=field*n+cells;
  yp=y;yp(cols)=yp(cols)+h(cols);[fp,gp]=plt.rhs(t,yp,c,p);
  df=fp-f;dg=gp-g;
  for j=1:numel(cells)
   cell=cells(j);col=cols(j);near=max(1,cell-1):min(n,cell+1);
   rows=reshape(bsxfun(@plus,near(:),(0:4)*n),[],1);
   values=[df(rows);dg(:,cell)]/h(col);rows=[rows;(5*n+1:m)'];
   ii=cursor+(1:numel(rows));ri(ii)=rows;ci(ii)=col;vv(ii)=values;
   cursor=cursor+numel(rows);
  end
 end
end
% Keep a fixed stored pattern for sparse solver backends at uniform states.
vv(vv==0)=1e-30;
J=sparse(ri(1:cursor),ci(1:cursor),vv(1:cursor),m,m);
end

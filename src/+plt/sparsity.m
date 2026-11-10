function J = sparsity(n)
%SPARSITY Nearest-neighbour five-field stencil plus diagnostic integrals.
B=spdiags(ones(n,3),-1:1,n,n);J=sparse(5*n+7,5*n+7);
J(1:5*n,1:5*n)=kron(ones(5),B);J(5*n+1:end,1:5*n)=1;
end

function P = pcg_basis(A,b,maxit,B)
%% primal formualtion
N = size(A,1);
uold = zeros(N,1);
rold = b - A*uold;
zold = B(rold);
pold = zold;    

P = zeros(N,maxit+1);
P(:,1) = pold/sqrt(pold'*A*pold);

for iter=1:maxit
    Ap = A*pold;

    alphaold = (rold'*zold)/(pold'*Ap);
    unew = uold + alphaold*pold;
    rnew = rold - alphaold*Ap;
    znew = B(rnew);

    beta = (rnew'*znew)/(rold'*zold);
    pnew = znew + beta*pold;
    P(:,iter+1) = pnew/sqrt(pnew'*A*pnew);
    
    uold = unew;
    pold = pnew;
    rold = rnew;
    zold = znew;
end







%% Solve the fractional equation (-\Delta)^s u = f by reduced basis  

% best uniform rational approximation
%res = [0.0002689 0.0055848 0.0272036 0.0965749 0.3202068 2.5105702];
%pol = [0.0 0.0000122 0.0006621 0.0127955 0.1626313 3.2129222];

res = [0.001163181622200 0.013218620357233 0.000255953763580 0.046695455564245 0.002833475749545 0.000144842262876 7.333904150551077 0.000535823283956 0.176401964217136 0.000085661432787 0.026987066006317 0.005823374267035 0.418860226358461 0.000717641159456 0.001870543346405 0.063983537529953 0.000318817071900 0.000126652196126 -0.002953457313354 0.494928246538564];
pol = [0.000018062500000 0.000933302500000 0.000000562500000 0.013156090000000 0.000104040000000 0.000000062500000 25.000000000000000 0.000003422500000 0.136641122500000 0.000000002500000 0.003271840000000 0.000277222500000 0.684011702500000 0.000008122500000 0.000044222500000 0.039243610000000 0.000001440000000 0.000000202500000 0.001759802500000 2.185962250000000];
np = length(res);

maxit = 20; fvH = cell(1); Pro = cell(maxit,1);
[node, elem] = cubemesh([0, 1, 0, 1, 0, 1],1);
bdFlag = [1 0 0 1;1 0 0 1;1 0 0 1;1 0 0 1;1 0 0 1;1 0 0 1];
[node,elem,bdFlag,fvH{1}] = mguniformrefine3(node,elem,bdFlag);
for iter=1:5
    Nold=size(node,1);
    [node,elem,bdFlag,fvH{iter+1},HB] = mguniformrefine3(node,elem,bdFlag);
    Nnew=size(node,1);
    
    Pro{iter}=sparse(1:Nold,1:Nold,1,Nnew,Nold);
    Pro{iter}=Pro{iter}+sparse(HB(:,1),HB(:,2),0.5,Nnew,Nold);
    Pro{iter}=Pro{iter}+sparse(HB(:,1),HB(:,3),0.5,Nnew,Nold);
end

NT = size(elem,1);
NV = size(node,1);
fprintf("The problem size is %d\n", NV);

[D, vol] = gradbasis3(node, elem);
[A, M] = P1mat3d(D,vol,elem);
rhs = P1rhs3d(node, elem, vol, @f);

T = auxmesh3(elem);
bv = unique(T.bdFace);
fv = setdiff(1:NV, bv);

bF = rhs(fv);
AF = A(fv, fv);
MF = M(fv, fv);
stiff_rb = AF + MF;

lambda = 18/vol(1).^(2/3);
bF = bF/lambda^0.5;
 
[Ai,Pro] = mgsetup(AF,Pro,fvH,100,1);
Bfun = @(x0) mgalg(x0,Ai,Pro,'V');

tic;
num_basis = 10;
P = pcg_basis(stiff_rb,rhs(fv),num_basis,Bfun);
t = toc;
fprintf('Elapsed time for reduced bases generation is %f seconds\n', t);

erri = zeros(num_basis,np);
ruhi = zeros(NV,np);
for i=1:np
    stiff_i = AF/lambda + pol(i)*MF;
    uhi = zeros(NV,1);
    [uhi(fv),cgerr,iter] = mypcg(stiff_i,bF,1e-6,1000,Bfun);
    fprintf('The number of iterations is %d, the error is %e\n',iter,cgerr);
    for m=1:num_basis
        Puh = ( P(:,1:m)'*stiff_i*P(:,1:m) )\( P(:,1:m)'*bF );
        ruhi(fv,i) = P(:,1:m)*Puh;
        erri(m,i) = sqrt( (uhi-ruhi(:,i))'*M*(uhi-ruhi(:,i)) );
        %fprintf('The rbm with %d bases has error %e\n', m, sqrt( (uh-ruh)'*(A+M)*(uh-ruh) ) );
    end
end
ruh = zeros(NV,1);
for i=1:np
    ruh = ruh + res(i)*ruhi(:,i);
end
l2 = CheckError3(ruh, node, elem, vol, @u);
uh = u(node); 
h1 = sqrt( (uh-ruh)'*A*(uh-ruh) );
fprintf('The rbm with %d bases has error %e\t%e\n', num_basis, l2, h1);
%--------------------------------------------------------------------------
function z = f(p)
% load data (right hand side function)
z = (3*pi^2)^0.5*u(p);
end
%--------------------------------------------------------------------------
function z = u(p) 
% exact solution of the test problem
z = sin(pi*p(:,1)).*sin(pi*p(:,2)).*sin(pi*p(:,3));
end
%--------------------------------------------------------------------------
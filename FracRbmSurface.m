%% Solve the fractional equation (-\Delta)^s u = f by reduced basis on 2d surfaces 
res = [0.001163181622200 0.013218620357233 0.000255953763580 0.046695455564245 0.002833475749545 0.000144842262876 7.333904150551077 0.000535823283956 0.176401964217136 0.000085661432787 0.026987066006317 0.005823374267035 0.418860226358461 0.000717641159456 0.001870543346405 0.063983537529953 0.000318817071900 0.000126652196126 -0.002953457313354 0.494928246538564];
pol = [0.000018062500000 0.000933302500000 0.000000562500000 0.013156090000000 0.000104040000000 0.000000062500000 25.000000000000000 0.000003422500000 0.136641122500000 0.000000002500000 0.003271840000000 0.000277222500000 0.684011702500000 0.000008122500000 0.000044222500000 0.039243610000000 0.000001440000000 0.000000202500000 0.001759802500000 2.185962250000000];
np = length(res);

node = [1, 0, 0; 0, 1, 0; -1, 0, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1];
elem = [1, 2, 5;3 5 2; 3, 4, 5; 1, 5, 4;...
        1, 6, 2;3 2 6; 3, 6, 4; 1, 4, 6];
for i=1:8
    [node, elem] = uniformrefine(node, elem);
end
node = node./sqrt(sum(node.*node,2));

NT = size(elem,1);
NV = size(node,1);
fprintf("The problem size is %d\n", NV);

[D, area] = gradbasis_surface(node, elem);
[A, M] = P1mat2d(D,area,elem);
rhs = P1rhs2d(node, elem, area, @f);

T = auxstructure(elem);

stiff_rb = A + 2*M;

lambda = 14/area(1);
Alambda = A/lambda;
rhs = rhs/lambda^0.5;

[Ai,Pro] = amgsetup( A + M );
Bfun = @(x0) mgalg(x0,Ai,Pro,'W');

tic;
num_basis = 9;
P = pcg_basis(stiff_rb,rhs,num_basis,Bfun);
t = toc;
fprintf('Elapsed time for reduced bases generation is %f seconds\n', t);

err = zeros(np,1);
ruhi = zeros(NV,np);
uhi = ruhi;
for i=1:np
    stiff_i = Alambda + pol(i)*M;
    uhi(:,i) = stiff_i\rhs;
    Puh = ( P'*stiff_i*P )\( P'*rhs );
    ruhi(:,i) = P*Puh;
    err(i) = sqrt( (uhi(:,i)-ruhi(:,i))'*M*(uhi(:,i)-ruhi(:,i)) );
    fprintf('The rbm for the %d-th pole has error %e\n', i, err(i) );
end
uh = u(node);
ruh = zeros(NV,1);
for i=1:np
    ruh = ruh + res(i)*ruhi(:,i);
end
fprintf('\n The rbm with %d bases has error %e\n', num_basis+1, sqrt( (uh-ruh)'*M*(uh-ruh) ) );
%showsolution(node,elem,ruh);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Sub functions called by main
%--------------------------------------------------------------------------
function z = f(p)
% load data (right hand side function)
z = (2)^0.5*u(p);
end
%--------------------------------------------------------------------------
function z = u(p) 
% exact solution of the test problem
z = p(:,1)+2*p(:,2)+3*p(:,3);
end
%--------------------------------------------------------------------------

%% approximate the function z^s on [1e-6,1] with iter OGA iterations
function [res, pol, err] = OGA(s,iter)
f = @(z) z.^(-s);
a = 1e-8; b = 1; 
node = unique([a:(0.001-a)/2000:0.001,0.001:(0.01-0.001)/2000:0.01,0.01:(b-0.01)/3000:b]');
h = node(2:end) - node(1:end-1); % graded mesh on [a,b]

c = [1/2-sqrt(15)/10 1/2 1/2+sqrt(15)/10];% 3-point Gauss quadrature
qpt = [node(1:end-1)+c(1)*h;node(1:end-1)+c(2)*h;node(1:end-1)+c(3)*h];

hd = 5e-5; t = (hd:hd:5)'; t = t.^2; % discretization of the dictionary
nd = length(t);
normg = sqrt(1./(a+t)-1./(b+t)); % L2 norm of dictionary elements
g = 1./(repmat(qpt,1,nd)+t')./normg'; % normalizing dictionary elements
  
fqpt = f(qpt); r = fqpt;
A = zeros(iter,iter); rhs = zeros(iter,1); % matrix and rhs for projection
id = zeros(iter,1); argmax = zeros(nd,1);
for i = 1:iter
    for j = 1:nd
        argmax(j) = product(g(:,j),r,h);
    end
    [~,id(i)] = max(abs(argmax));
    for j = 1:i
        A(j,i) = product(g(:,id(j)),g(:,id(i)),h);
        A(i,j) = A(j,i);
    end
    rhs(i) = product(g(:,id(i)),fqpt,h);
    C = lsqminnorm(A(1:i,1:i),rhs(1:i));
    r = fqpt - g(:,id(1:i))*C;
    fprintf("Step %d\n",i);
end
res = C./normg(id); % residues of the approximant 
pol = t(id); % (-1)*poles of the approximant
maxnode = linspace(1e-6,1,5000000)';
rf = zeros(length(maxnode),1);
for j=1:iter
    rf = rf + res(j)./(maxnode+pol(j));
end
err = max(abs(rf - f(maxnode))); % max norm error on [1e-6,1]
%--------------------------------------------------------------------------
function z = product(f1,f2,h)
z = 5/18*sum( h.*f1(1:3:end).*f2(1:3:end) )+...
    4/9*sum( h.*f1(2:3:end).*f2(2:3:end) )+...
    5/18*sum( h.*f1(3:3:end).*f2(3:3:end) );
end

end



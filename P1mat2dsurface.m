function [A,M] = P1mat2dsurface(D,area,elem)
%% Assemble P1 matrices A, M  
NV = max(elem(:));
A = zeros(NV,NV);
for i = 1:3
    for j = 1:3
        Aij = area.*dot(D(:,:,i).*D(:,:,j),2);
        A = A + sparse(elem(:,i), elem(:,j), Aij, NV, NV);
    end
end
ii=[elem(:,1) elem(:,2) elem(:,3) elem(:,1) elem(:,1) elem(:,2) elem(:,2) elem(:,3) elem(:,3)];
jj=[elem(:,1) elem(:,2) elem(:,3) elem(:,2) elem(:,3) elem(:,3) elem(:,1) elem(:,1) elem(:,2)];
M = sparse(ii(:),jj(:),[repmat(area/6,3,1);repmat(area/12,6,1)],NV,NV);




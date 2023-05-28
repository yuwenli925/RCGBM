function [A,M] = P1mat2d(D,area,elem)
%% Assemble P1 matrices A, M  
NV = max(elem(:));

ii=[elem(:,1) elem(:,2) elem(:,3) elem(:,1) elem(:,1) elem(:,2) elem(:,2) elem(:,3) elem(:,3)];
jj=[elem(:,1) elem(:,2) elem(:,3) elem(:,2) elem(:,3) elem(:,3) elem(:,1) elem(:,1) elem(:,2)];
M = sparse(ii(:),jj(:),[repmat(area/6,3,1);repmat(area/12,6,1)],NV,NV);

D12 = area.*dot(D(:,:,1),D(:,:,2),2);
D13 = area.*dot(D(:,:,1),D(:,:,3),2);
D23 = area.*dot(D(:,:,2),D(:,:,3),2);
A = sparse(ii(:),jj(:),[area.*dot(D(:,:,1),D(:,:,1),2);area.*dot(D(:,:,2),D(:,:,2),2);area.*dot(D(:,:,3),D(:,:,3),2);...
    D12;D13;D23;D12;D13;D23],NV,NV);



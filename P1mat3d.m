function [A, M] = P1mat3d(D,vol,elem)
%% assemble the P1 stiffness matrix in 3d
NV = max(elem(:));

A=sparse(elem(:),elem(:),[vol.*dot(D(:,:,1),D(:,:,1),2);...
    vol.*dot(D(:,:,2),D(:,:,2),2);vol.*dot(D(:,:,3),D(:,:,3),2);...
    vol.*dot(D(:,:,4),D(:,:,4),2)],NV,NV);
M=sparse(elem(:),elem(:),repmat(vol/10,4,1),NV,NV);

ii=[elem(:,1);elem(:,1);elem(:,1);elem(:,2);elem(:,2);elem(:,3)];
jj=[elem(:,2);elem(:,3);elem(:,4);elem(:,3);elem(:,4);elem(:,4)];
U = sparse(ii,jj,repmat(vol/20,6,1),NV,NV);
M = M + U + U';

U = sparse(ii,jj,[vol.*dot(D(:,:,1),D(:,:,2),2);vol.*dot(D(:,:,1),D(:,:,3),2)
    vol.*dot(D(:,:,1),D(:,:,4),2);vol.*dot(D(:,:,2),D(:,:,3),2);...
    vol.*dot(D(:,:,2),D(:,:,4),2);vol.*dot(D(:,:,3),D(:,:,4),2)],NV,NV);
A = A + U + U';

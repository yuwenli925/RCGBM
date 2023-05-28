function [D,area,normal] = gradbasis_surface(node,elem)
%% GRADBASIS gradient of barycentric basis. 

ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);

normal = cross(ve3,ve1,2); area2 = vecnorm(normal,2,2);
area = area2/2;

normal = normal./(area*2);

D(:,:,1) = cross(normal,ve1,2)./area2;
D(:,:,2) = cross(normal,ve2,2)./area2;
D(:,:,3) = cross(normal,ve3,2)./area2;

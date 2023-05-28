function rhs = P1rhs2d(node, elem, area, f)

fpxy = f( ( node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:) )/3 );   
bt = repmat(fpxy.*area/3,1,3);
rhs = accumarray(elem(:),bt(:),[size(node,1) 1]);
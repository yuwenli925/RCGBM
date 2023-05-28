function rhs = P1rhs3d(node, elem, vol, f)

fpxy = f( ( node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:)+node(elem(:,4),:) )/4 );   
bt = repmat(fpxy.*vol/4,1,4);
rhs = accumarray(elem(:),bt(:),[size(node,1) 1]);
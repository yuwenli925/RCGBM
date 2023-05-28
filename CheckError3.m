function l2 = CheckError3(uh, node, elem, vol, u)

[lambda, weight] = quadpts3(2);
num = length(weight);
l2err = zeros(size(elem,1),1); 
h1err = l2err;
for p=1:num

    uhpxy = uh(elem(:,1))*lambda(p,1) ...
          + uh(elem(:,2))*lambda(p,2) ...
          + uh(elem(:,3))*lambda(p,3) ...
          + uh(elem(:,4))*lambda(p,4);

    point = node(elem(:,1),:)*lambda(p,1) ...
          + node(elem(:,2),:)*lambda(p,2) ...
          + node(elem(:,3),:)*lambda(p,3) ...
          + node(elem(:,4),:)*lambda(p,4);

    l2err = l2err + weight(p)*( u(point) - uhpxy ).^2;
end   

l2 = sqrt(vol'*l2err);

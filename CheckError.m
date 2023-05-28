function error = CheckError(uh, node, elem, vol, u)

[lambda, weight] = quadpts(3);
num = length(weight);
l2error = zeros(size(elem,1),1); 

for p=1:num

    uhpxy = uh(elem(:,1))*lambda(p,1) ...
          + uh(elem(:,2))*lambda(p,2) ...
          + uh(elem(:,3))*lambda(p,3);

    point = node(elem(:,1),:)*lambda(p,1) ...
          + node(elem(:,2),:)*lambda(p,2) ...
          + node(elem(:,3),:)*lambda(p,3);

    l2error = l2error + weight(p)*( u(point) - uhpxy ).^2;
end   

error = sqrt(vol'*l2error);

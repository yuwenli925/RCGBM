function [node,elem,bdFlag,fv,HB] = mguniformrefine3(node,elem,bdFlag)
%% UNIFORMREFINE3 uniformly refine a 3-D triangulation.
% 
% [node,elem] = uniformrefine3(node,elem) divides each tetrahedra into
% eight small similar sub-tetrahedrons.
%
% [node,elem,bdFlag,fv] = uniformrefine3(node,elem,bdFlag) also update boundary
% conditions represented by bdFlag.  fv is the set of free vertices.
%
% [node,elem,bdFlag,fv,HB] = uniformrefine3(node,elem) outputs HB array 
% which is useful for nodal interpolation.

if ~exist('bdFlag','var'), bdFlag =[]; fv = []; end
if ~exist('HB','var'),    HB = [];     end

%% Construct data structure
[elem2dof,edge] = dof3P2(elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Add new nodes
node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2;   
HB(:,[1 2 3]) = [(N+1:N+NE)', edge(:,1:2)]; 


%% Refine each tetrahedron into 8 tetrahedrons
t = 1:NT;
p(t,1:10) = elem2dof;
elem(8*NT,:) = [0 0 0 0]; % enlarge the elem array
elem(t,:) = [p(t,1), p(t,5), p(t,6), p(t,7)];
elem(NT+1:2*NT,:) = [p(t,5), p(t,2), p(t,8), p(t,9)];
elem(2*NT+1:3*NT,:) = [p(t,6), p(t,8), p(t,3), p(t,10)];
elem(3*NT+1:4*NT,:) = [p(t,7), p(t,9), p(t,10), p(t,4)];
% always use diagonal 6-9. The ordering is important. See the reference.
elem(4*NT+1:5*NT,:) = [p(t,5), p(t,6), p(t,7), p(t,9)];
elem(5*NT+1:6*NT,:) = [p(t,5), p(t,6), p(t,8), p(t,9)];
elem(6*NT+1:7*NT,:) = [p(t,6), p(t,7), p(t,9), p(t,10)];
elem(7*NT+1:8*NT,:) = [p(t,6), p(t,8), p(t,9), p(t,10)];

%% Update boundary edges
if ~isempty(bdFlag)
    bdFlag(8*NT,:) = [0 0 0 0]; % enlarge the bdFlag array
    bdFlag(NT+1:2*NT,[1 3 4]) = bdFlag(t,[1 3 4]); 
    bdFlag(2*NT+1:3*NT,[1 2 4]) = bdFlag(t,[1 2 4]); 
    bdFlag(3*NT+1:4*NT,[1 2 3]) = bdFlag(t,[1 2 3]);
    % always use diagonal 6-9
    bdFlag(4*NT+1:5*NT,2) = bdFlag(t,3);
    bdFlag(5*NT+1:6*NT,4) = bdFlag(t,4);
    bdFlag(6*NT+1:7*NT,3) = bdFlag(t,2);
    bdFlag(7*NT+1:8*NT,1) = bdFlag(t,1);
    % change t in the last
    bdFlag(t,1) = 0;
    
    bdEI = sum(bdFlag,2); bdEI = ~(bdEI==0);
    bdElem = elem(bdEI,:); 
    bdF = bdFlag(bdEI,:); 
    bdvec = ~bdF(:); el2v = bdElem(:);
    bv = el2v(bdvec); NV = max(elem(:));
    fv = setdiff((1:NV)',bv);
else
    bdFlag = [];
end

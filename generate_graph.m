function [L,Adj,Weights,Div]=generate_graph(n)

%% given n, this generates a graph Laplacian with random weights; 
%% It also generates the discrete gradient:

L=triu(sprand(n,n,5/n),1); %% 5/n is the density; 
[i,j,s]=find(L+L');
s=rand(length(s),1);
L=sparse(i,j,s,n,n)+spdiags([ones(n,1),ones(n,1)],[-1,1],n,n);
L=0.5*(L+L');
[i,j,s]=find(L);
o=ones(length(s),1);
Adj=sparse(i,j,o,n,n); %%adjacency matrix:
%% make now $L$ a graph laplacian:
L=spdiags([sum(L)'],[0],n,n)-L;
%% the grad maps space of edges to the space of vertices
[iv,jv,Weights]=find(triu(L,1));
nedges=length(iv);
Weights=-spdiags(Weights,[0],nedges,nedges);

Div=sparse(iv,[1:nedges]',ones(nedges,1),n,nedges);
Div=Div+sparse(jv,[1:nedges]',-ones(nedges,1),n,nedges);
%%Check if all is OK:
% chk=norm(Div*Weights*Div'-L,inf)
%% A good indefinite problem would be with a matrix L1, defined as
%% L1=[D0,Div';Div,sparse(n,n)]; %% here D0 is an (nedges by nedges SPD matrix could be diagonal). 
end
function [Ai,Pro] = mgsetup(A1,Pro1,fv1,N0,dirichlet)
%% Construct transfer operator
Pro1=Pro1(~cellfun('isempty',Pro1));
level = length(Pro1)+1; Pro=Pro1;
if (level==1)||(size(A1,1)<N0)
    Ai{1}=A1;Pro=cell(1);
else
    %level = floor(level/m)+1;
    Ai = cell(level,1);
    Ai{level}=A1;
    for k=level-1:-1:1
        if size(Ai{k+1},1)>N0
            if dirichlet
                Pro{k}=Pro1{k}(fv1{k+1},fv1{k});
                Ai{k} = Pro{k}'*Ai{k+1}*Pro{k};
            else
                Pro{k}=Pro1{k};
                Ai{k} = Pro{k}'*Ai{k+1}*Pro{k};
            end
        else
            Ai{k}=[];Pro{k}=[];
        end
    end
    Ai=Ai(~cellfun('isempty',Ai));
    Pro=Pro(~cellfun('isempty',Pro));
    Pro{end+1}=[];
end




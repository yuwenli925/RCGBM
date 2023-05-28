function [Ai,Pro,cl] = amgsetup(A)
%% classical algebraic multigrid
interpolation = 'S'; theta = 0.025; N0 = 100;
%% Construct transfer operator
N = size(A,1);
level = max(min(ceil(log2(N)/2-4),8),2);
Ai = cell(level,1);
Pro = cell(level,1);
Res = cell(level,1);
Ai{level} = A;
cl = 1; % coarsest level
for k = level:-1:2
    % coarsening
    [isC,As] = coarsenAMGc(Ai{k},theta);
    % prolongation and restriction
    switch interpolation
        case 'S'  % standard
            [tempPro,tempRes] = interpolationAMGs(As,isC);
        case 'T'  % two points
            [tempPro,tempRes] = interpolationAMGt(As,isC);
    end
    % record operators   
    Pro{k-1} = tempPro;
    Res{k} = tempRes;
    Ai{k-1} = tempRes*Ai{k}*tempPro;
    % check if reach the coarsest level
    if size(Ai{k-1},1)< N0 
        cl = k-1;
        Ai=Ai(~cellfun('isempty',Ai));
        Pro=Pro(~cellfun('isempty',Pro));
        Pro{end+1}=[];
        break;
    end
end

end



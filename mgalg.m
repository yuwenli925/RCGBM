function x = mgalg(b,Ai0,Pro0,solver)
%% geometric multigrid
mu = 1;
level = length(Ai0); Ai = Ai0; Pro = Pro0;
Nb = length(b); NA = size(Ai{end},1);
N = Nb/NA; x = zeros(Nb,1);

%% V-cycle and W-cycle
if solver=='V'
    for i=1:N
        ii=(i-1)*NA+1:i*NA;
        x(ii) = vcycle(b(ii),level);
    end
end

if solver=='W'
    x = wcycle(b,level);
end

if strcmp(solver,'BPX')
    x = bpx(b,level);
end
%% V cycle recursive function
function e = vcycle(r,J)        % solve equations Ae = r in each level
    if J == 1
        e = Ai{1}\r;   % exact solver in the coaresest grid
        return
    end
    % fine grid pre-smoothing
    e = triu(Ai{J})\r;   % pre-smoothing
    for s = 1:mu           % extra mu steps smoothing
        e = e + triu(Ai{J})\(r-Ai{J}*e); 
    end
    rc = Pro{J-1}'*(r - Ai{J}*e);
    % coarse grid correction twice
    ec = vcycle(rc,J-1);
    % fine grid post-smoothing
    e = e + Pro{J-1}*ec;
    e = e + tril(Ai{J})\(r-Ai{J}*e);
    for s = 1:mu
        e = e + tril(Ai{J})\(r-Ai{J}*e); % post-smoothing
    end
end

%% W cycle recursive function
function e = wcycle(r,J)        % solve equations Ae = r in each level
    if J == 1
        e = Ai{1}\r;   % exact solver in the coaresest grid
        return
    end
    % fine grid pre-smoothing
    e = triu(Ai{J})\r;   % pre-smoothing
    for s = 1:mu           % extra mu steps smoothing
        e = e + triu(Ai{J})\(r-Ai{J}*e); 
    end
    rc = Pro{J-1}'*(r - Ai{J}*e);
    % coarse grid correction twice
    ec = wcycle(rc,J-1);
    ec = ec + wcycle(rc - Ai{J-1}*ec,J-1);
    % fine grid post-smoothing
    e = e + Pro{J-1}*ec;
    e = e + tril(Ai{J})\(r-Ai{J}*e);
    for s = 1:mu
        e = e + tril(Ai{J})\(r-Ai{J}*e); % post-smoothing
    end
end

%% BPX 
function e = bpx(r,J)        % solve equations Ae = r in each level
    if J == 1
        %e = r./diag(Ai{1});   % exact solver in the coaresest grid
        e = Ai{1}\r;
        return
    end
    % fine grid pre-smoothing
    e = r./diag(Ai{J});
    % coarse grid correction 
    ec = bpx(Pro{J-1}'*r,J-1);
    e = e + Pro{J-1}*ec;
end
end


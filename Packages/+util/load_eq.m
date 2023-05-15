function state = load_eq(domain,eq_name)

    % loads eq from saved source
%     if(eq_name(1)~='R')
%         eq_name = strcat('R40_',eq_name);
%     end
%     data = load(strcat('Equilibria/',eq_name,'.mat'));
    data = load([eq_name,'.mat']);
    if(isfield(data,'s'))
        state = data.s;
    else
        u = data.u1;
        v = data.u2;
        state(:,:,1) = u;
        state(:,:,2) = v;
    end
    
    % resizes eq to domain size
    [n,m,~] = size(state);
    p = domain.Nx;
    q = domain.Ny;
    state = util.fftgridchange(state,n,m,p,q);
end


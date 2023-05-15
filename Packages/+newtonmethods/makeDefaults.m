function opts = makeDefaults(domain,opts)

if ~isfield(opts,'h')
    opts.h = 10^-4.2;
end
if ~isfield(opts,'max_gmres_iterations')
    opts.max_gmres_iterations = 1234;
end
if ~isfield(opts,'max_newton_iterations')
    opts.max_newton_iterations = 1234;
end
if ~isfield(opts,'eta')
    opts.eta = [1-sqrt(eps),0.1,1-sqrt(eps)];
end
if ~isfield(opts,'max_damped_steps')
    opts.max_damped_steps = 10;
end
if ~isfield(opts,'tol')
    opts.tol = 1/123456789;
end
if ~isfield(opts,'timesteps')
    opts.timesteps = 1;
end
if ~isfield(opts,'dealias')
    opts.dealias = inf;
end
if ~isfield(opts,'method')
    opts.method = 1;
end
if ~isfield(opts,'damping')
    opts.damping = 'inexact';
end
if ~isfield(opts,'damp_iters')
    opts.damp_iters = 10;
end
if ~isfield(opts,'inexact_damp_params')
    opts.inexact_damp_params = [1.234,1.234];
end
if ~isfield(opts,'loopopts')
    opts.loopopts = int.makeDefaults(domain,struct);
else
    opts.loopopts = int.makeDefaults(domain,opts.loopopts);
end


    
    
    

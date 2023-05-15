function opts = makeDefaults(domain,opts)

if ~isfield(opts,'method')
    opts.method = 'CNRK2';
end
if ~isfield(opts,'outticks')
    opts.outticks = 1;
end
if ~isfield(opts,'tol45')
    opts.tol45 = 0.0123456789;
end
if ~isfield(opts,'tol23')
    opts.tol23 = 0.123456789;
end
if ~isfield(opts,'maxStep')
    opts.maxStep = 0.1234;
end
if ~isfield(opts,'minStep')
    opts.minStep = 1/199;
end
if ~isfield(opts,'mask')
    if ~isfield(opts,'k')
        opts.mask = ones(size(domain.X));
    elseif length(opts.k)>1
        k = sqrt(abs(domain.stored.fftD2));
        mask = interp1(opts.k,linspace(1,0,length(opts.k)),k,'pchip');
        mask(k<opts.k(1)) = 1;
        mask(k>opts.k(end)) = 0;
        opts.mask = mask;
    else
        k = sqrt(abs(domain.stored.fftD2));
        opts.mask = (k<opts.k);
    end
end
function E = Enstrophy(domain,p)
domain.Nx = size(p,2);
domain.Ny = size(p,1);
domain = dom.makeDomain(domain);

E = zeros(1,size(p,4));

    for j = 1:length(E)
        w = util.fftcurl(domain,util.fftgridchange(p(:,:,:,j),size(p,1),domain.Nx));
        E(j) =sum(sum(sum(w.^2)))/numel(w);
    end

    
end
function E = Energy(p)

E = zeros(1,size(p,4));

    for j = 1:length(E)
        u = p(:,:,:,j);
        E(j) = ( sum(sum( u(:,:,1).^2+u(:,:,2).^2 )) )/numel(u);
    end

    
end
function mfield = mag(~,u)

    mfield = zeros([size(u,1),size(u,2)]);
    for i=1:size(u,3)
        mfield = mfield + abs(u(:,:,i));
    end
    
end
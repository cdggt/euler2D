function vis_1dfield(domain,field)
warning('off', 'MATLAB:contour:ConstantData')
% domain -- domain object that the data was created on
% field  -- can either be a single velocity field or a trajectory itself

x = domain.X;
y = domain.Y;
xbnd = [x(1) x(end)];
ybnd = [y(1) y(end)];

figure 

h = contourf(x,y,field, 20,'EdgeColor','None');
colormap parula
colorbar
xlim(xbnd);
ylim(ybnd);

warning('on', 'MATLAB:contour:ConstantData')
end


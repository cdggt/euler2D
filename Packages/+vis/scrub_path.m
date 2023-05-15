function [] = scrub_path(domain,path)

A = 1;
B = size(path,4);

uf = uifigure;

ax = axes('Parent',uf,'position',[0.1 0.3  0.8 0.6]);
[~,h] = contourf(ax,util.fftstream(domain,path(:,:,:,1)),30,'Edgecolor','none');
colormap(ax,'parula');

b = uislider(uf);

b.Limits = [A B];
b.Position = [20 80 500 3];          
b.ValueChangingFcn = {@callback,h,domain,path};
b.MajorTicks = linspace(A,B,B-A+1);

end

function [] = callback(es,ed,handle,domain,path)
    handle.ZData = util.fftstream(domain,path(:,:,:,round(ed.Value)));
end
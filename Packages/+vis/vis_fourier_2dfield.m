function vis_fourier_2dfield(~,field,title_label_x,title_label_y)

figure 

U = fftshift(fft2(field(:,:,1)));
V = fftshift(fft2(field(:,:,2)));
subplot(2,1,1)
imagesc(abs(U));
if(exist('title_label_x','var'))
    title(title_label_x,'Interpreter','Latex');
end
colorbar
subplot(2,1,2)
imagesc(abs(V));
if(exist('title_label_y','var'))
    title(title_label_y,'Interpreter','Latex');
end
colorbar

dim = get(gcf,'position');
set(gcf,'position',[dim(1),dim(2),dim(3),2*dim(4)]);

end
function vis_fourier_1dfield(~,field,title_label)

figure 

Field = fftshift(fft2(field));

imagesc(abs(Field));
if(exist('title_label','var'))
    title(title_label,'Interpreter','Latex');
end
colorbar

end
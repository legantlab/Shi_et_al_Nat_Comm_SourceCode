function img_GT = gen_multibead3D_GT(img_size,wavelength,ri,density)
pixel_size = wavelength*0.1/ri/1000;         %pixel size in um, PSF are based on 0.1 media wavelength per pixel
edge_pix =5;
img_area = (img_size(1)-2*edge_pix)*(img_size(2)-2*edge_pix)*(img_size(3)-2*edge_pix)*(pixel_size^3);
n_bead = round(img_area*density);
bead_rad = 0.05/pixel_size;
% bead_rad = 0.2/pixel_size;

% x_range = (img_size(2)-1)/2;
% z_range = (img_size(1)-1)/2;


[XX,YY,ZZ] = meshgrid(1:1:img_size(2),1:1:img_size(1),1:1:img_size(3));

coord = [edge_pix edge_pix edge_pix] + rand(n_bead,3).*[(img_size(2)-2*edge_pix) , (img_size(1)-2*edge_pix), (img_size(3)-2*edge_pix)];

img_GT = zeros(img_size);

for i =1:n_bead
    img_GT(sqrt((XX-coord(i,1)).^2+(YY-coord(i,2)).^2+(ZZ-coord(i,3)).^2)<=bead_rad) =  img_GT(sqrt((XX-coord(i,1)).^2+(YY-coord(i,2)).^2+(ZZ-coord(i,3)).^2)<=bead_rad) + 1;
end


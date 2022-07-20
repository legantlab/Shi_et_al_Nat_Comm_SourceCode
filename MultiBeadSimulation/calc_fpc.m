function [flc_out,flc_area] = calc_fpc(img1,img2)
%Calculate Fourier plane correlation based on source codes from Nieuwenhuizen et al(2013) 
%Input: img1, img2: 3D image stacks of two independent images of taken on
%the same field of view.
%Output: flc_out: Fourier plane correlation of the pair images. Y axis is the axial direction and X axis is the lateral axis. The output assume
%homogenous resolution in the lateral plane.
%flc_area: integrated area with value > 1/7 in Fourier plane correlation

% Removing edge artifacts in Fourier space using Periodic plus Smooth Image Decomposition
[img1_p,~] = perdecomp_3D(img1);
[img2_p,~] = perdecomp_3D(img2);

img1_dip = dip_image(img1_p);
img2_dip = dip_image(img2_p);

raw_window = round(size(img1,1)/2.4);
center_pt = (size(img1,1)+1)/2;


[flc_out_dip_xy,flc_out_dip_xz,flc_out_dip_yz]= fpc(img1_dip,img2_dip);
flc_out_dip = (flc_out_dip_xz+flc_out_dip_yz)/2;
flc_out = dip_array(flc_out_dip);
% flc_out_sub  = flc_out(center_pt-raw_window:center_pt+raw_window,center_pt-raw_window:center_pt+raw_window);
% flc_area = sum(flc_out_sub>1/7,'all');
stats = regionprops(flc_out>1/7,'centroid','Area');
flc_area = max([stats.Area]);
end


function output_PSF = filter_PSF_3d(ref_PSF,input_PSF)
%FILTER_PSF Summary of this function goes here
%   Detailed explanation goes here
n_img = size(input_PSF,1);
threshold = 1E-3;
for i =1 :n_img
    ref_OTF = fftshift(fftn(ifftshift(squeeze(ref_PSF(i,:,:,:)))));
    input_OTF = fftshift(fftn(ifftshift(squeeze(input_PSF(i,:,:,:)))));
    input_OTF(abs(ref_OTF)<max(abs(ref_OTF(:)))*threshold)=0;
    output_PSF(i,:,:,:) = abs(fftshift(ifftn(ifftshift(input_OTF))));

end


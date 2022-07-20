function img_conv = conv3d_fast(img,psf_input)
pad_size = 200;
otf_input = fftshift(fftn(ifftshift(psf_input)));
img_fft = fftshift(fftn(ifftshift(img)));
img_conv = abs(fftshift(ifftn(ifftshift(otf_input.*img_fft))));
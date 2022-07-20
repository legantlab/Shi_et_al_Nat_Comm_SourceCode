%-- needs the matlab toolbox dipimage, free download at www.diplib.org
clear;
root_folder = '..\'; % folder of the source code


% path for prerequsites
addpath('..\FIREfunctions');   % folder for FIRE files
run('C:\Program Files\DIPimage 2.9\dipstart.m');                  % run DIP, follow instructions from www.diplib.org
det_PSF_path = 'Z:\Yu\PSFGenerator_detPSF\detPSF\RW_det_NA1p0_RI_1p3_WL560_PSF.mat';    % path to 3D detection PSF
PSF_folder = 'Z:\Yu\2021_05_26(LLS_beam_characterization_paper)\SourceCodeForNatComm\BeamSimulation\simulation_output\Test\Gaussian_newNA_p21p0_fill_factor_1_propagation_80';   % folder for PSF 


print_flag = true;


bead_intensity = 1E3;     % Intensity of the bead on average
wavelength = 560;   % wavelength in nm
ri = 1.3;       % refrective index
density =3;   % bead density in beads/um^2
img_size = [201 201 201];  % size of simulated images, in X, Y and Z
SBR = 10;      % SBR ratio of the image

n_deconiter = 10;   % Number of RL decon iterations

%% Load detection PSF 
if ~exist('psf','var')
    load(det_PSF_path);
end

%% Generate Overall PSF in 3D
output_folder = [PSF_folder '\multibead3D'];
mkdir(output_folder);
% load in excitation PSF
PSF_filename = [PSF_folder '\LLS_data.mat'];
load(PSF_filename,'NormalizedDitheredxzPSF','DitheredIntensityz0');

% Calculate the positon of FWHM along the propagation direction
Normalized_intensity = DitheredIntensityz0/DitheredIntensityz0(1);
[~,fwhm_index] = min(abs(Normalized_intensity-0.5));
fwhm_halfindex = round((1+fwhm_index)/2);

OverallDitheredPSF = [];

% Calculate Overall PSF as a product of detection and excitation PSF
for i_plane = [1 fwhm_halfindex fwhm_index]
    CurrentePSF = squeeze(NormalizedDitheredxzPSF(i_plane,:,:));
    CurrentePSF_3D = repmat(CurrentePSF,[1 1 size(NormalizedDitheredxzPSF,2)]);
    CurrentePSF_3D = permute(CurrentePSF_3D,[3,2,1]);
    OverallDitheredPSF= cat(4,OverallDitheredPSF,CurrentePSF_3D.*psf);
end
OverallDitheredPSF = permute(OverallDitheredPSF,[4 1 2 3]);
PSF_mid_plane = (size(OverallDitheredPSF,2)+1)/2;
half_size = (img_size(1)-1)/2;


OverallDitheredPSF_crop = OverallDitheredPSF(:,PSF_mid_plane-half_size:PSF_mid_plane+half_size,PSF_mid_plane-half_size:PSF_mid_plane+half_size,PSF_mid_plane-half_size:PSF_mid_plane+half_size);

PSF_size = size(OverallDitheredPSF_crop);

% Normalize PSF to have a peak intensity of 1
OverallDitheredPSF_normalized = OverallDitheredPSF_crop./sum(OverallDitheredPSF_crop,[2:4]);

%% Generate ground truth image
GT_folder = [root_folder '\multibead_groudTruth'];
mkdir(GT_folder);


GT_name = [replace(['3d_wl' num2str(wavelength) '_ri' num2str(ri) '_density' num2str(density) '_size' num2str(img_size(1)) '_' num2str(img_size(1))],'.','p') '.mat'];
GT_file = dir([GT_folder '\' GT_name]);

if isempty(GT_file)
    bead_GT = gen_multibead3D_GT(img_size,wavelength,ri,density);
    save([GT_folder '\' GT_name],'bead_GT');
else
    load([GT_file(1).folder '\' GT_file(1).name]);
end


amplifi_ratio_PSF = 3E7;
bead_signal = bead_intensity*300;
bead_GT = bead_GT*bead_signal;

noise_level = bead_signal/600/SBR;    % Background level in ground truth image
bead_GT1 = bead_GT + 2*noise_level;%+noise_level*randn(size(bead_GT));
bead_GT1(bead_GT1<0)=0;
bead_GT2 = bead_GT + 2*noise_level;%+noise_level*randn(size(bead_GT));
bead_GT2(bead_GT2<0)=0;

%% Simulate Raw multibead images using loaded PSF
file_name = ['Multibead Density' num2str(density,'%1.0e') ' SNR' num2str(SBR,'%1.0e')];


% simulate 2 instances of raw multibead image at center with Poisson noise
multibead_center1 = poissrnd(conv3d_fast(bead_GT1,squeeze(OverallDitheredPSF_normalized(1,:,:,:))));
multibead_center1(multibead_center1<0)=0;
multibead_center2 = poissrnd(conv3d_fast(bead_GT2,squeeze(OverallDitheredPSF_normalized(1,:,:,:))));
multibead_center2(multibead_center2<0)=0;
% calculate FPC based on the two simulated images 
[multibead_flc_center,multibead_flcarea_center] = calc_fpc(multibead_center1,multibead_center2);

% simulate 2 instances of raw multibead image at half way between center and FWHM with Poisson noise
multibead_sidehalf1 = poissrnd(conv3d_fast(bead_GT1,squeeze(OverallDitheredPSF_normalized(2,:,:,:))));
multibead_sidehalf1(multibead_sidehalf1<0)=0;
multibead_sidehalf2 = poissrnd(conv3d_fast(bead_GT2,squeeze(OverallDitheredPSF_normalized(2,:,:,:))));
multibead_sidehalf2(multibead_sidehalf2<0)=0;
% calculate FPC based on the two simulated images 
[multibead_flc_sidehalf,multibead_flcarea_sidehalf] = calc_fpc(multibead_sidehalf1,multibead_sidehalf2);

% simulate 2 instances of raw multibead image at FWHM with Poisson noise
multibead_side1 = poissrnd(conv3d_fast(bead_GT1,squeeze(OverallDitheredPSF_normalized(3,:,:,:))));
multibead_side1(multibead_side1<0)=0;
multibead_side2 = poissrnd(conv3d_fast(bead_GT2,squeeze(OverallDitheredPSF_normalized(3,:,:,:))));
multibead_side2(multibead_side2<0)=0;
% calculate FPC based on the two simulated images 
[multibead_flc_side,multibead_flcarea_side] = calc_fpc(multibead_side1,multibead_side2);
%RMSE_side = sqrt(mean((multibead_side(:)-bead_GT(:)).^2));


x_ticks = [1:(img_size(2)-1)/4:img_size(2)];
x_label = [-(img_size(2)-1)/2:(img_size(2)-1)/4:(img_size(2)-1)/2]*0.1;
z_ticks = [1:(img_size(1)-1)/4:img_size(1)];
z_label = [-(img_size(1)-1)/2:(img_size(1)-1)/4:(img_size(1)-1)/2]*0.1;


flc_winsize = (img_size(1)-1)/4;
x_midpt = (img_size(2)+1)/2;
z_midpt = (img_size(1)+1)/2;
flcx_ticks = [1:(img_size(1)-1)/8:1+2*flc_winsize];
flcx_label = {'-1/4','-1/8','0','1/8','1/4'};
flcz_ticks = [1:(img_size(2)-1)/8:1+2*flc_winsize];
flcz_label = {'-1/4','-1/8','0','1/8','1/4'};





figure('Position',[0 0 1200 300]);
subplot(1,4,1)
imagesc(imrotate(squeeze(bead_GT((img_size(1)+1)/2,:,:)),90));
axis equal
title('Ground Truth')
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');

subplot(1,4,2)
imagesc(imrotate(squeeze(multibead_center1((img_size(1)+1)/2,:,:)),90));
axis equal
title(['image at focus, FRC area:' num2str(multibead_flcarea_center)])
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');

subplot(1,4,3)
imagesc(imrotate(squeeze(multibead_sidehalf1((img_size(1)+1)/2,:,:)),90));
axis equal
title(['image at halfway FWHM, FRC area:' num2str(multibead_flcarea_sidehalf)])
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');

subplot(1,4,4)
imagesc(imrotate(squeeze(multibead_side1((img_size(1)+1)/2,:,:)),90));
axis equal
title(['image at FWHM, FRC area:' num2str(multibead_flcarea_side)])    
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');

sgtitle(file_name);
if print_flag
    saveas(gcf,[output_folder '\' file_name '.tif']);
end

figure('Position',[0 0 800 200]);
subplot(1,3,1)
imagesc(multibead_flc_center(z_midpt-flc_winsize:z_midpt+flc_winsize,x_midpt-flc_winsize:x_midpt+flc_winsize))
axis equal
caxis([0 1]);
title(['FLC at focus, FRC area:' num2str(multibead_flcarea_center)])
xticks(flcx_ticks);
xticklabels(flcx_label);
xlabel('kx (1/\lambda)');
yticks(flcz_ticks);
yticklabels(flcz_label);
ylabel('kz (1/\lambda)');


subplot(1,3,2)
imagesc(multibead_flc_sidehalf(z_midpt-flc_winsize:z_midpt+flc_winsize,x_midpt-flc_winsize:x_midpt+flc_winsize))
axis equal
caxis([0 1]);
title(['FLC at halfway FWHM, FRC area:' num2str(multibead_flcarea_sidehalf)])
xticks(flcx_ticks);
xticklabels(flcx_label);
xlabel('kx (1/\lambda)');
yticks(flcz_ticks);
yticklabels(flcz_label);
ylabel('kz (1/\lambda)');

subplot(1,3,3)
imagesc(multibead_flc_side(z_midpt-flc_winsize:z_midpt+flc_winsize,x_midpt-flc_winsize:x_midpt+flc_winsize))
axis equal
caxis([0 1]);
title(['FLC at FWHM, FRC area:' num2str(multibead_flcarea_side)])  
xticks(flcx_ticks);
xticklabels(flcx_label);
xlabel('kx (1/\lambda)');
yticks(flcz_ticks);
yticklabels(flcz_label);
ylabel('kz (1/\lambda)');

sgtitle(file_name);

if print_flag 
    saveas(gcf,[output_folder '\flc_' file_name '.tif']);
end

flc_area_raw = [multibead_flcarea_center,multibead_flcarea_sidehalf,multibead_flcarea_side];

%% RL decon

PSF_noise = 1E-6*amplifi_ratio_PSF;

multibead_flc_center_RLdecon_mean = zeros(size(bead_GT,2:3));
multibead_flc_sidehalf_RLdecon_mean = zeros(size(bead_GT,2:3));
multibead_flc_side_RLdecon_mean = zeros(size(bead_GT,2:3));

% add noise to simulated PSF 
OverallDitheredPSF_normalized_1p = poissrnd(OverallDitheredPSF_normalized*amplifi_ratio_PSF) +PSF_noise*randn(size(OverallDitheredPSF_normalized));
OverallDitheredPSF_normalized_1p(OverallDitheredPSF_normalized_1p<0)=0;
% Filter out components in frequency space out of OTF support
OverallDitheredPSF_normalized_1p = filter_PSF_3d(OverallDitheredPSF_crop,OverallDitheredPSF_normalized_1p);
OverallDitheredPSF_normalized_1p(OverallDitheredPSF_normalized_1p<0)=0;
OverallDitheredPSF_normalized_1p = OverallDitheredPSF_normalized_1p./sum(OverallDitheredPSF_normalized_1p,[2:4]);




flcx_ticks_decon = [1:(img_size(2)-1)/4:img_size(2)];
flcx_label_decon = {'-1/2','-1/4','0','1/4','1/2'};
flcz_ticks_decon = [1:(img_size(1)-1)/4:img_size(1)];
flcz_label_decon = {'-1/2','-1/4','0','1/4','1/2'};


multibead_center_RLdecon1 = deconvlucy(multibead_center1,squeeze(OverallDitheredPSF_normalized_1p(1,:,:,:)),n_deconiter);
multibead_center_RLdecon1(multibead_center_RLdecon1<0)=0;
multibead_sidehalf_RLdecon1 = deconvlucy(multibead_sidehalf1,squeeze(OverallDitheredPSF_normalized_1p(2,:,:,:)),n_deconiter);
multibead_sidehalf_RLdecon1(multibead_sidehalf_RLdecon1<0)=0;
multibead_side_RLdecon1 = deconvlucy(multibead_side1,squeeze(OverallDitheredPSF_normalized_1p(3,:,:,:)),n_deconiter);
multibead_side_RLdecon1(multibead_side_RLdecon1<0)=0;


figure('Position',[0 0 1200 300]);
subplot(1,4,1)
imagesc(imrotate(squeeze(bead_GT((img_size(1)+1)/2,:,:)),90));
axis equal
title('Ground Truth')
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');


subplot(1,4,2)
imagesc(imrotate(squeeze(multibead_center_RLdecon1((img_size(1)+1)/2,:,:)),90));
axis equal
title(['image at focus'])
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');

subplot(1,4,3)
imagesc(imrotate(squeeze(multibead_sidehalf_RLdecon1((img_size(1)+1)/2,:,:)),90));
axis equal
title(['image at halfway FWHM'])
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');

subplot(1,4,4)
imagesc(imrotate(squeeze(multibead_side_RLdecon1((img_size(1)+1)/2,:,:)),90));
axis equal
title(['image at FWHM'])  
xticks(x_ticks);
xticklabels(x_label);
xlabel('x (\lambda)');
yticks(z_ticks);
yticklabels(z_label);
ylabel('z (\lambda)');

sgtitle([file_name ' decon RL']);
if print_flag
    saveas(gcf,[output_folder '\deconRL_' file_name '.tif']);
end

save([output_folder '\' file_name '.mat']);
    
close all;

save([output_folder '\summary.mat'],'flc_area_raw');

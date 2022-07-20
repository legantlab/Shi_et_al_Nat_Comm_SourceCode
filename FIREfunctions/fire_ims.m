%FIRE_IMS   Compute resolution from two images
%
% The FRC curve and subsequently the resolution are computed from the two
% input images.
%
% SYNOPSIS:
%   [fire_value frc_out fire_high fire_low] = fire_ims(in1,in2,pixelsize,show_frc)
%
%   pixelsize
%      Pixel size of the images (in nm)
%   show_frc
%      Display FRC curve?
%
% NOTES:
%   Non-square images are zero padded to make them square.
%
% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen & Bernd Rieger, Dec 2012

function varargout = fire_ims(varargin)

d = struct('menu','FIRE',...
           'display','FIRE from images',...
           'inparams',struct('name',       {'in1',      'in2',          'SR_pixelsize',     'show_frc'},...
                             'description',{'Image 1',  'Image 2',      'Pixel size (nm)',  'Display FRC curve'},...
                             'type',       {'image',    'image',        'array',            'boolean'},...
                             'dim_check',  {2,          2,              0,                  0},...
                             'range_check',{[],         [],             [eps Inf],          []},...
                             'required',   {1,          1,              1,                  0},...
                             'default',    {'in1',      'in2',          10,                  1}...
                              ),...
           'outparams',struct('name',{'fire_value','frc_curve','fire_high','fire_low'},...
                              'description',{'FIRE value','FRC curve','FIRE - 1 std. dev.','FIRE + 1 std. dev.'},...
                              'type',{'array','array','array','array'},...
                              'suppress',{1,1,1,1}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
      varargout{1} = d;
      return
   end
end

try
   [in1, in2, SR_pixelsize, show_frc] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end

%% Compute results
[fire_value, frc_curve, fire_high, fire_low] = frcresolution(in1,in2);

fprintf('FIRE value %2.1f +- %2.2f nm.\n', fire_value*SR_pixelsize, (fire_low-fire_high)/2*SR_pixelsize);
fprintf('FIRE value %2.1f +- %2.2f pixels.\n', fire_value, (fire_low-fire_high)/2);

%% Plot FRC curve
if show_frc
    qmax = 0.5/(SR_pixelsize);
    
    figure
    hold on
    plot([0 qmax],[0 0],'k')
    plot(linspace(0,qmax*sqrt(2), length(frc_curve)), frc_curve,'-')
    plot([0 qmax],[1/7 1/7],'m')
    plot(1/(fire_value*SR_pixelsize),1/7,'rx')
    plot(1/(fire_value*SR_pixelsize)*[1 1],[-0.2 1/7],'r')
    hold off
    xlim([0,qmax]);
    ylim([-0.2 1.2])
    xlabel('Spatial frequency (nm^{-1})');
    ylabel('FRC')
end

%% Outputs
varargout{1} = fire_value*SR_pixelsize;
varargout{2} = frc_curve;
varargout{3} = fire_high*SR_pixelsize;
varargout{4} = fire_low*SR_pixelsize;
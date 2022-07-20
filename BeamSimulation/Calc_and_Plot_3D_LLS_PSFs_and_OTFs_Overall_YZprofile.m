function [propagation_length, DitheredIntensityz0, PWb] = Calc_and_Plot_3D_LLS_PSFs_and_OTFs_Overall(lattice_descrip, xyPol, PW, NAmax, NAideal, NAmin, NAdet, fill_factor, crop_factor, numysteps, ystepsize, xz_det_PSF, xz_det_OTF, gamma,root_folder,namin_ratio);
%
%Calc_and_Plot_3D_LLS_PSFs_and_OTFs 
%   The lattice is calculated using the algorithm where the ideal 2D 
%   lattice is bound in z, and the real part of the bound
%   lattice E field is written to the SLM.  This creates an E field
%   impinging on the mask that has equal kz bounding for all k vectors.
%   This field is then filtered by an annular mask of NAmax ans NAmin
%   determined by the input values of these variables, and approximately
%   equal to the kz bounding when fill_factor = 1.

%
%   After determining the post mask E field, the program then calculates
%   the overall xz PSFs and OTFs at various positions along the propagation
%   direction y to see how the PSFs and the relative amplitudes of the 
%   excitation spatial frequencies vary form their initital amplitudes on focus
%
%   propagation_length is an estimate of the +/- extent of the light sheet 
%      in the y direction, and should roughy equal the 50% intensity points
%      of the sheet
%   DitheredIntensityz0 is a numysteps length vector giving the
%      intensity of the lattice at the z = 0 central plane at various positions 
%      y along the propagation vector
%   PWb is an N x 7 matrix where the first six columns are PW from above,
%      and the 7th column is the fill factor ratios between the different k
%      vectors
%
%IF lattice_descrip = "Gaussian"
%   a.  PW is overwritten with PW = [1 0 0 0 0 0] to let the program know
%       the only k vector is one along the axis of the objective
%   b.  A circular mask of NA = 0.6 is used, with no DC block (since the Gaussian beam is DC) 
%END
%
%   PW is an N x 6 matrix containing the directions of the N plane
%      waves that constitute the lattice.  Each row refers to a single 
%      plane wave.  The properties of each plane wave are entered in the
%      columns of PW as follows:
%      PW(n, 1:3) = complex electric field vector.  
%         IMPORTANT:  FOR THIS PROGRAM ONLY, PW(n,1) contains the relative
%         E field amplitude that should be used for each k vector, and 
%         PW(n,2:3) are ignored, since the final E field compoments
%         PW(n,1:3) are calculated within the program itself
%      PW(n, 4:6) = direction ek
%   NAmax is the maximum numerical aperture in the rear pupil of the
%       excitation objective
%   NAideal is the numerical aperture of the ideal 2D lattice
%       NAmax >= NAideal >= NAmin
%   NAmin is the minimum numerical aperture in the rear pupil of the
%       excitation objective
%   fill_factor is a positive number giving the desired ratio between the
%      half-width of the illumination at the Bessel annulus and the width of
%      the annulus itself, for the case where all k vectors have the same
%      bounding
%   numysteps is the number of steps along the propagation axis y at which
%      the LLS should be characterized
%   ystepsize is the y distance between these characterization steps along
%      y in media wavelengths
%
if strcmp(lattice_descrip, 'Gaussian')
   PW = [1 0 0 0 0 0];  %enforce a single k vector at the center of the objective for Gaussian
   NAmax = NAmax - NAmin
   NAideal = 0;
   NAmin = 0;   %enforces a centrally located Gaussian lone of pupil 1/e^2
                %half-width = NAmax
end
%
index = 1.33;
namax = NAmax./index;
%For MB square lattice, the NA of the two side beamlets will be different
%from the two at kx = 0
if strcmp(lattice_descrip,'Bessel')
    naideal = zeros(size(PW,1),1);
    naideal(3:4) = NAideal./index;
    naideal(1:2) = NAmin./index.*namin_ratio;
else
    naideal = NAideal./index;
end
namin = NAmin./index;
%
%calc the rear pupil E field for this nominal condition where all k vectors 
%   have the same bounding:
NominalPupilEfield = zeros(1001,1001); %match the pupil field of other 
%   programs that use pixsize = 0.1 media wavelengthsand +/-50 media 
%   wavelgnths field of view. In this case, the pupil field calculation 
%   runs from +/-2.5 times the max spatial frequency defined by two
%   counterpropagating plane waves
pixsize = 0.1;
plot_range = 50;
%
%estimate the lattice propagation length, based on the largest value of 
%  max-min illumination NA among the various bounded k vectors of the lattice:
%first calculate the nominal kz gaussian envelope of the Efield by having it
%   fall to 1/e^2 in amplitude at the edges of the annulus, for the case of
%   fill_factor = 1.  Then, for all fill_factors:
if strcmp(lattice_descrip,'Bessel')
    kzsigma = (namax - namin).*fill_factor;  %E = Eo*exp(-(kz/(2*kzsigma))^2)
    kzsigma_side = (sqrt(namax^2 - naideal(1)^2)).*fill_factor;
else
    kzsigma = (namax - namin).*fill_factor;  %E = Eo*exp(-(kz/(2*kzsigma))^2)
end
%now find the range dNA of illumination associated with the +/-kz edges
%   of each bound k vector:
B = size(PW);
dNA = zeros(1,B(1));
namax_beamlet = dNA;
namin_beamlet = dNA;
if strcmp(lattice_descrip,'Bessel')
    for p = 1:B(1)
        namax_beamlet(p) = sqrt(naideal(p).^2 + 2.*naideal(p).*kzsigma.*abs(PW(p,5)) + kzsigma.^2);    
        if naideal(p).*abs(PW(p,5)) >= kzsigma  
            %the current illumination stripe does not pass the equator 
            namin_beamlet(p) = sqrt(naideal(p).^2 - 2.*naideal(p).*kzsigma.*abs(PW(p,5)) + kzsigma.^2);
        else
            %the current illumination stripe does intersect the equator
            namax_beamlet(p) = sqrt(naideal(p).^2 + 2.*naideal(p).*kzsigma.*abs(PW(p,5)) + kzsigma.^2);
            namin_beamlet(p) = naideal(p);
        end
        dNA(p) = namax_beamlet(p) - namin_beamlet(p);
    end
else
    for p = 1:B(1)
        namax_beamlet(p) = sqrt(naideal.^2 + 2.*naideal.*kzsigma.*abs(PW(p,5)) + kzsigma.^2);
        if naideal.*abs(PW(p,5)) >= kzsigma  
            %the current illumination stripe does not pass the equator
            namin_beamlet(p) = sqrt(naideal.^2 - 2.*naideal.*kzsigma.*abs(PW(p,5)) + kzsigma.^2);
        else
            %the current illumination stripe does intersect the equator
            namin_beamlet(p) = naideal;
        end
        dNA(p) = namax_beamlet(p) - namin_beamlet(p);
    end    
end

dNA
%
%find the minimum dNA among the k vectors, corresponding to the beam with
%   the least bounding and hence largest beam waist at the focus:
[dNAmin, kindex] = min(dNA);

if strcmp(lattice_descrip,'Bessel')
    kymax = 2.*pi.*sqrt(1-namin.^2);
    kymin = 2.*pi.*sqrt(1-namax.^2);    
else
    kymax = 2.*pi.*sqrt(1-namin_beamlet(kindex).^2);
    kymin = 2.*pi.*sqrt(1-namax_beamlet(kindex).^2); 
end

kydiff = kymax - kymin;
propagation_length = 2.*pi./kydiff  %approximate full +/-extent of the lattice in y, in media wavelengths

folder = [root_folder '\' lattice_descrip '_NA_p' num2str(round(100*NAmax)) 'p' num2str(round(100*NAmin)) '_fill_factor_' num2str(fill_factor) '_propagation_' num2str(round(propagation_length)) '\'];
mkdir(folder);
%
%plot at multiple points along the propagation axis:
yp = ystepsize.*(0:1:(numysteps-1));
%
%find the period in x and z of the ideal lattice:
%   each period is given by the smallest non-zero difference of the
%      k vector components along that axis
%   max. confinement in each direction is given by the largest
%      difference of the k vector components along that axis
B = size(PW);
for m = 1:2   %loop thru both axes
    mindiff = 2;   %difference must be < 2 for normalized k
    maxdiff = 0;
    for n = 1:B(1)
        for q = 1:B(1)  %double loop thru all plane wave pairs
            PWdiff = abs(PW(n,(m + 3)) - PW(q,(m + 3)));  % abs of the sum
                        %of the mth component of the nth and qth k vectors
            if ((PWdiff > 0.001) & (PWdiff < mindiff))
                %if difference non-zero yet smaller than the smallest
                %difference thus far then
                mindiff = PWdiff;
            end
            if PWdiff > maxdiff
                %if difference is larger than the largest
                %difference thus far then
                maxdiff = PWdiff;
            end
        end
    end   %done looping thru all plane wave pairs
    period(1,m) = 1 ./ mindiff;  %period in wavelengths along m axis
    period(2,m) = 1 ./ maxdiff;  %confinement in wavelengths along m axis
end
%
if strcmp(lattice_descrip,'Bessel')
    %here making kzsigma to fill the entire na gap, both beamlets at
    %equator and at center, with equal intensity throughout the whole beamlet
    kzsigma = (namax - namin).*fill_factor/2;  %E = Eo*exp(-(kz/(2*kzsigma))^2)
    kzsigma_side = (sqrt(namax^2 - naideal(1)^2)).*fill_factor;    
    %
    B = size(PW);
    for p = 1:B(1)
      %find the ideal center location of each k vector:
      centerx = 501 + round(250.*naideal(p).*PW(p,4)./2.5);
      centerz = 501 + round(250.*naideal(p).*PW(p,5)./2.5);
      EAmp = PW(p,1);  %manual adjustment of k vector stored in first index of PW
    %express kzsigma in terms of pixels:
      if p < 3
          kzsigma_pix = 250.*kzsigma_side./2.5;
      else
          kzsigma_pix = 250.*kzsigma./2.5;
      end
      EFieldkz = EAmp.*double(abs((1:1:1001)-centerz)<kzsigma_pix)';                  % change to a flat top function instead of Gaussian
      NominalPupilEfield(:,centerx) = NominalPupilEfield(:,centerx) + EFieldkz;
    end
    NominalPupilIntensity = NominalPupilEfield.*conj(NominalPupilEfield);
    NominalPupilIntensity = NominalPupilIntensity./max(max(NominalPupilIntensity));    
    %
    %plot the nominal intensity in the rear pupil before adjustment:
    NominalPupilIntensity2 = NominalPupilIntensity(300:700,300:700);
    C = size(NominalPupilIntensity2);    
    %
    %     %calc and store the ratios of dNA for each k vector, normalized to the
    %     %dNAmin minimum value among all k vectors, as the 7th index in an 
    %     %expanded PW matrix:
    PWb = zeros(B(1),7);
    PWb(:,1:6) = PW;
    PWb(:,7) = dNAmin./dNA;

    AdjustedPupilEfield = NominalPupilEfield;
    ESLM = fftshift(ifft2(ifftshift(AdjustedPupilEfield)));
    RealE = real(ESLM);
    maxval = max(max(RealE));
    minval = min(min(RealE));
    PlotRealE = (RealE - minval) ./ (maxval - minval);
    A = size(PlotRealE);
    %
    %calculate and plot the SLM pattern that hopefully gives rise to the
    %   desired E field at the pupil:
    RealE = RealE ./ (max(max(abs(RealE))));
    SLMPhase = RealE'.*pi;
    
else
    %insure that the k vectors lie in a plane (i.e., cone angle = 90 deg):
    PW(:,4:5) = PW(:,4:5) ./ sqrt(1 - PW(1,6)^2);
    PW(:,6) = 0;
    %
    %now modify each k vector to reflect the cone angle upon which they lie:
    PW(:,4:5) = naideal .* PW(:,4:5);
    PW(:,6) = sqrt(1-naideal.^2);
    %
    %normalize the input polarization and make it a 1 x 3 vector:
    InputPol = [xyPol(1) xyPol(2) 0];
    InputPol = InputPol ./ norm(InputPol);
    %
    %find the electric field for each lattice wavevector after passing through
    %the objective:
    B = size(PW);
    for n = 1:B(1)
        %find the orthonormal vectors defining the coordinate system for the
        %nth beam when entering the rear pupil:
        phivec = cross(squeeze(PW(n,4:6)), [0 0 1]);
        phivec = phivec ./ norm(phivec);  %azimuthal unit vector
        radvec = cross([0 0 1], phivec);  %radial unit vector
        %
        %the azimuthal component of the electric field is unaffected when passing
        %through the objective:
        ephi = dot(phivec, InputPol);
        %
        %the radial component is tilted by refraction when passing through the
        %objective to point in the theta direction as defined by a spherical
        %coordinate system centered at the focal point:
        thetavec = cross(squeeze(PW(n,4:6)), phivec);
        etheta = dot(radvec, InputPol);
        %
        %save the desired electric field amplitude stored in PW(n,1) before 
        %loading the complexelectric field into PW(n,1:3)
        EAmp = PW(n,1);
        %the sum of the azimuthal and theta components gives the total electric
        %field for the nth plane wave of the lattice:
        PW(n,1:3) = ephi .* phivec + etheta .* thetavec;
        %
        %confirm that the electric field is of unit strength:
        PW(n,1:3) = PW(n,1:3) ./ norm(PW(n,1:3));
        %adjust the eletric field to the desired amplitude:
        PW(n,1:3) = EAmp.*PW(n,1:3);
    end
    %
    %calculate the complete electric field of the ideal 2D lattice:
    x = 0:pixsize:plot_range;
    y = 0:pixsize:plot_range;
    [X Y] = ndgrid(x, y);
    %now calculate the E field at each mesh point:
    A = size(X);
    Ex = zeros(A(1), A(2));
    Ey = zeros(A(1), A(2));
    Ez = zeros(A(1), A(2));
    for q = 1:B(1)   %loop thru all plane waves
        phase = exp(2 .* pi .* i .* (PW(q, 4) .* X + PW(q, 5) .* Y));
        Ex = Ex + PW(q, 1) .* phase;
        Ey = Ey + PW(q, 2) .* phase;
        Ez = Ez + PW(q, 3) .* phase;
    end
    %expand through all quadrants:
    Extemp = zeros(2 .* A(1), 2 .* A(2));
    Eytemp = zeros(2 .* A(1), 2 .* A(2));
    Eztemp = zeros(2 .* A(1), 2 .* A(2));
    %load the original data into the first quadrant:
    Extemp((A(1) + 1):end, (A(2) + 1):end) = Ex;
    Eytemp((A(1) + 1):end, (A(2) + 1):end) = Ey;
    Eztemp((A(1) + 1):end, (A(2) + 1):end) = Ez;
    %now mirror along each dimension and use parities to fill other quadrants:
    %simply mirror the data since parity is always even for magnitudes:
    Extemp(1:A(1), (A(2) + 1):end) = flipdim(Ex, 1);
    Eytemp(1:A(1), (A(2) + 1):end) = flipdim(Ey, 1);
    Eztemp(1:A(1), (A(2) + 1):end) = flipdim(Ez, 1);
    Extemp(:, 1:A(2)) = flipdim(Extemp(:, (A(2) + 1):end), 2);
    Eytemp(:, 1:A(2)) = flipdim(Eytemp(:, (A(2) + 1):end), 2);
    Eztemp(:, 1:A(2)) = flipdim(Eztemp(:, (A(2) + 1):end), 2);
    %delete the extra vector from mirroring in each dimension:
    Extemp(A(1), :) = [];
    Eytemp(A(1), :) = [];
    Eztemp(A(1), :) = [];
    Extemp(:, A(2)) = [];
    Eytemp(:, A(2)) = [];
    Eztemp(:, A(2)) = [];
    Ex = Extemp;
    Ey = Eytemp;
    Ez = Eztemp;
    clear Extemp Eytemp Eztemp;
    A = size(Ex);
    %
    %find the ideal 2D lattice intensity:  %%%%%%
    EComp = zeros(3,A(1),A(2));
    EComp(1, :, :) = Ex;
    EComp(2, :, :) = Ey;
    EComp(3, :, :) = Ez;
    ESq = EComp .* conj(EComp);
    ESqTot = squeeze(sum(ESq, 1));
    maxval = max(max(ESqTot));
    ESqTot = ESqTot ./ maxval;
    %
    %calc and plot the ideal 2D lattice of the component of the real electric field 
    %projected onto the state of the input electric field polarization:
    RealE = real(conj(Ex) .* InputPol(1) + conj(Ey) .* InputPol(2));
    %                                  
    %calc and plot the bound 2D electric field lattice at the SLM:
    z = -plot_range:pixsize:plot_range;
    kxmax = 2.*pi.*namax;
    kxmin = 2.*pi.*namin;
    kxdiff = kxmax - kxmin;
    lattice_full_width = pi./kxdiff;  %approximate half width of the function limiting
        %the extent of the bound lattice, in media wavelengths
    sigma = lattice_full_width ./ sqrt(2.*log(2)) ./ fill_factor;
    envelope = exp(-2.*(z./sigma).^2)' * ones(1,A(1));
    RealE = RealE .* envelope';
    maxval = max(max(RealE));
    minval = min(min(RealE));
    PlotRealE = (RealE - minval) ./ (maxval - minval);
    %
    %now calc and plot the binary phase pattern at SLM needed to create the bound lattice
    %   at the sample:
    %first truncate any values below the level of crop_factor:
    RealE = RealE ./ (max(max(RealE)));
    RealE = RealE .* (abs(RealE) > crop_factor);
    %now create the binary phase function (0 or pi) for the SLM:
    SLMPhase = (sign(RealE + eps).*pi./2 + pi./2);
    PWb = PW;
end
%

%
%now find the pupil field at the mask created by this SLM pattern:
if strcmp(lattice_descrip, 'Gaussian')
    kzsigma = (namax - namin).*fill_factor;  %E = Eo*exp(-(kz/(2*kzsigma))^2)
    %express kzsigma in terms of pixels:
    kzsigma_pix = 250.*kzsigma./2.5;
    EImpingingOnMask = zeros(1001,1001);
    EImpingingOnMask(501,:) = exp(-((-500:1:500)./kzsigma_pix).^2);
    A=size(EImpingingOnMask);
else
    maxval = max(max(SLMPhase));
    minval = min(min(SLMPhase));
    PlotSLMPhase = (SLMPhase' - minval) ./ (maxval - minval);
    
    ESLM = exp(1i .* SLMPhase);
    EImpingingOnMask = fftshift(fft2(ifftshift(ESLM)));  %complex electric field impinging on the annular mask
    %there is a singularity in EImpiningOnMask at center pixel 501,501 which we remove by
    %   replacing it with the value at the adjacent pixel 501,500:
    EImpingingOnMask(501,501) = EImpingingOnMask(501,500);
end
%plot the intensity impinging on the annular mask:
IntensityOnMask = EImpingingOnMask .* conj(EImpingingOnMask);   %intensity impinging on the annular mask
maxval = max(max(IntensityOnMask));
IntensityOnMask = IntensityOnMask ./ maxval;
%crop the result to the +/- the max spatial freq defined by
%conterpropagating plane waves:
IntensityOnMask = IntensityOnMask(300:700,300:700);
%
%calc and plot the function for the filtering provided by the annular mask:
halfpix = (A(1)-1)./2;
x = -halfpix:halfpix;
x = 5.*x./halfpix;   %match the spectral range of the electric field at the annular mask
y = x;
[X Y] = ndgrid(x, y);
R = sqrt(X.*X + Y.*Y);

%the pupil field is filtered by a thin annulus of size NAmin to NAmax
MaxRad = namax;  %maximum annulus diameter, normalized to index
MinRad = namin;  %minimum annulus diameter, normalized to index

if strcmp(lattice_descrip, 'Gaussian')
   MaxRad = 0.6./index;
   MinRad = 0;
elseif max(abs(PW(:,5))) == 0
   %if the k vectors are confined to y = 0, then use namax, namid for the
   %   annulus to block the bigger DC beam that exists in this case
   MaxRad = namax;  %maximum annulus diameter, normalized to index
   MinRad = namin;  %minimum annulus diameter, normalized to index   
end
AnnularFilter = (R <= MaxRad) & (R >= MinRad);
%
%calc the E field immediately after transmission through the annulus:
EAfterMask = EImpingingOnMask .* AnnularFilter;
%shift EAfterMask by one pixel to correctly center it within the pupil:
EAfterMask(:,1:1000) = EAfterMask(:,2:1001);
%plot the intensity immediately after the annular mask:
MaskIntensity = EAfterMask .* conj(EAfterMask);
MaskIntensity = MaskIntensity ./ max(max(MaskIntensity));
MaskIntensity = MaskIntensity(300:700,300:700);

%
%calc and plot the LLS xz excitation PSF at the focal plane of the excitation objective:
ESample = fftshift(ifft2(ifftshift(EAfterMask)));
SampleIntensity = ESample .* conj(ESample);
SampleIntensity = SampleIntensity ./ max(max(SampleIntensity));

%
%now calculate the LLS xz excitation pattern and overall OTF at points y<>0
%   along the propagation direction
%
%to speed the 2D integration over all points in the rear pupil, create a
%   list of all points in the pupil at which the amplitude of the
%   Efield is above some minimum threshold:
EAfterMask = EAfterMask';
EAmpAfterMask = abs(EAfterMask);
EAmpAfterMask = EAmpAfterMask./max(max(EAmpAfterMask));
EAmpThreshhold = 0.001;
% 
%loop through all rows and columns of EAmpAfterMask within the square
%   bounding the annulus:
halfrange = ceil(namax.*halfpix./2);
NumPointsWithinBoundingSquare = (2.*halfrange + 1).^2;
%define matrices to store the info about the illuminated points in the pupil:
NAx = zeros(1,NumPointsWithinBoundingSquare);
NAy = NAx;
PupilEfield = NAx;
numpoints = 1;
for p = (halfpix-halfrange):(halfpix+halfrange)
    for q = (halfpix-halfrange):(halfpix+halfrange)
        if EAmpAfterMask(p,q) > EAmpThreshhold
            NAx(numpoints) = p;
            NAy(numpoints) = q;
            PupilEfield(numpoints) = EAfterMask(p,q);
            numpoints = numpoints + 1;
        end
    end
end
numpoints = numpoints - 1
NAx = 2.5.*namax.*(NAx(1:numpoints) - halfpix)./halfrange;
NAy = 2.5.*namax.*(NAy(1:numpoints) - halfpix)./halfrange;
PupilEfield = PupilEfield(1:numpoints);
%
%find the k vector direction for all illuminated points in the pupil:
sinth = sqrt(NAx.^2 + NAy.^2);
costh = sqrt(1 - sinth.^2);  %angle of k vector to y axis
sinphi = NAy./sinth;
cosphi = NAx./sinth;
kx = cosphi.*sinth;
kz = sinphi.*sinth;
ky = costh;
%
%loop through all longitudinal points yp at which we wish to find the xz
%    excitation cross section of the light sheet as well as the xz overall
%    PSF and OTF:
ys = size(yp);
%
%define the plot parameters for the xz excitation at points y <> 0 along
%   the propagation direction:
xp = -50:0.1:50;
zp = xp;
%the above values are in media wavelengths
[XP, ZP] = ndgrid(xp, zp);
%
Q = size(XP);
XZIntensity = zeros(ys(2),Q(1),Q(2));  %stores the xz intensity at the 
                   %calculation points yp along the propagation direction
for p = 1:ys(2)
    EField = zeros(Q(1),Q(2));  %initialize the field before integration 
                                %over all illuminated points in the pupil
    for q = 1:numpoints
        EField = EField + PupilEfield(q).*exp(2 .* pi .* 1i .* (kx(q).*XP + kz(q).*ZP + ky(q).*yp(p)));
    end
    XZIntensity(p,:,:) = EField.*conj(EField);    
end
%normalize all the xz images along y to the maximum value, presumably at y = 0:
XZIntensity = XZIntensity./max(max(max(XZIntensity)));
%
%save YZ plane at x = 0
J = size(XZIntensity);
midzpix = (J(3)+1)./2;
SampleIntensity = imrotate(squeeze(XZIntensity(:,:,midzpix)),90);

figure('Color','white');  %create a new figure window for the plots
set(gcf, 'Position', [50 100 600 600]);
axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* SampleIntensity);
colormap hot(256);
axis([1 J(1) 1 J(2)]);
set(gca, 'XTick', [1:(J(1)-1)./5:J(1)]);
set(gca, 'XTickLabel', max(yp) .* [0:0.2:1]);
xlabel(['y / \lambda'], 'FontSize', 14);
set(gca, 'YTick', [1:(J(2)-1)./4:J(2)]);
set(gca, 'YTickLabel', max(zp) .* [-1:0.5:1]);
ylabel(['z / \lambda'], 'FontSize', 14);
text(0.05 .*J(2), -0.06 .* J(3), ['YZ Profile ', num2str(yp(p), '%1.2f'), '\lambda From Focus, NA ', num2str(NAmax, '%1.2f'), '/', num2str(NAideal, '%1.2f'), '/', num2str(NAmin, '%1.2f')], 'FontSize', 12);
text(0.15 .*J(2), -0.02 .* J(3), [lattice_descrip, ',  fill factor = ', num2str(fill_factor, '%1.2f')], 'FontSize', 12);
saveas(gcf,[folder 'yz_Profile_' num2str(round(yp(p))) 'wl.tif']);


save([folder 'LLS_data_yzprofile.mat'],'SampleIntensity','yp','zp'); 
end

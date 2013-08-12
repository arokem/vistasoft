function err = dtiRawRohdeEddyError(mc, phaseDir, srcIm, trgIm, sampDensity)
%  Minimization function for Rohode eddy/motion correction
% err = dtiRawRohdeEddyError(c, phaseDir, srcIm, trgIm, sampDensity)
%
% mc is a 1x14 array of the 6 motion params (translations, rotations) and
% the 8 eddy-correct params (c). 
%
% Apply the Rohde 14-parameter motion/eddy-current deformation to the
% source image and return normalized mutual information between the
% transformed source and target images.
%
% If the srcIm and trgIm  are not changing after the first call, you can pass
% an empty array for each of those since the joint-hist routine keeps a
% data cache. 
%
%   Rohde, Barnett, Basser, Marenco and Pierpaoli (2004). Comprehensive
%   Approach for Correction of Motion and Distortion in Diffusion-Weighted
%   MRI. MRM 51:103-114.
%
% TODO: The calculations should be done on physical-space coords (ie.
% scaled by mmPerVox).
%
%
% HISTORY:
%
% 2007.05.02 RFD wrote it.
% 2007.05.10 RFD: added data cache and in-place coord xform to
% mrAnatFastInterp3. That makes this function about 4-5x faster.

sz = size(trgIm);

mc = mc(:)';

% Allow rigid-body xform:
if(length(mc)==6) mc = [mc 0 0 0 0 0 0 0 0 0];
else mc = [mc phaseDir]; end

% Create the joint histogram
H = dtiJointHist(trgIm, srcIm, mc, sampDensity);

% Smooth the histogram
fwhm = 7;
lim  = ceil(2*fwhm);
s = (fwhm/sqrt(8*log(2)))^2+eps;
% Note- for small FWHM, the following will not be accurate.
krn = (1/sqrt(2*pi*s))*exp(-([-lim:lim].^2)/(2*s));
krn = krn./sum(krn);
H = conv2(krn,krn,H);

% Compute cost function
H  = H+eps;
sh = sum(H(:));
H  = H/sh;
s1 = sum(H,1);
s2 = sum(H,2);
% Normalised Mutual Information of:  Studholme,  Hill & Hawkes (1998).
% "A normalized entropy measure of 3-D medical image alignment".
% in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA,
% pp. 132-143.
err = (sum(s1.*log2(s1))+sum(s2.*log2(s2)))/sum(sum(H.*log2(H)));
% flip the sign so minimization maximizes NMI
err= -err;

return;

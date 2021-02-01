function [VolumeImage] = function_Rendering(Setup, z,DMDPatterns)
%This function renders images at distance z from DMD patterns
%Inputs:
%z = an array of size LZ corresponding to the propagation distance of the image depth planes
%DMDPatterns = a set of p DMD pattterns of size LX, by LY, by p
%Outputs:
%VolumeImage = A stack of intensity data corresponding to individually propagated oherent waves, of size 
%LX by LY  by LZ by p

[LX,LY,LP] = size(DMDPatterns);
UX = Setup.ps*(1:Setup.DMDX); UX = UX-mean(UX);
UY = Setup.ps*(1:Setup.DMDY); UY = UY-mean(UY);
[XX,YY] = ndgrid(UX,UY);
laser_amplitude = exp(-((XX.^2+YY.^2)/Setup.laserradius^2));
mask = double(sqrt(XX.^2+YY.^2)<0.0001);

LZ = numel(z);

VolumeImage = zeros(Setup.DMDX,Setup.DMDY,LZ,LP);
for k = 1:LP

FieldA = laser_amplitude.*squeeze(DMDPatterns(:,:,k));

%go to image plane
[FieldB,psx,psy] = function_lens(FieldA,Setup.ps,Setup.ps,-Setup.f,Setup.lambda);
 
FieldB = FieldB/sqrt((sum(abs(FieldB(:).^2))));

for j = 1:LZ

    fieldZ = function_propagate((1-mask).*FieldB,Setup.lambda,z(j),psy, psx);
   fieldZ = fieldZ/sqrt((sum(abs(fieldZ(:).^2))));
    VolumeImage(:,:,j,k) = abs(fieldZ.^2);
end
end

end

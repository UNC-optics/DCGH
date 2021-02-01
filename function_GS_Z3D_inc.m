function [DMDPattern] = function_GS_Z3D_inc(Setup, DMDPatternref,z,target_amplitude)
%This is a core step of the DCGH - Gerchberg–Saxton algorithm 
%to optimize several coherent wave modulation patterns at once.
%Inputs:
%target_amplitude = target 3D intensity distribution
%z = location of target along the z plane
%DMDPatternref = an intial estimation of the DMDPattern to be optimized 
%Outputs:
%DMDPattern = the updated DMD pattern


LP = numel(z);

UX = Setup.ps*(1:Setup.DMDX); UX = UX-mean(UX);
UY = Setup.ps*(1:Setup.DMDY); UY = UY-mean(UY);
[XX,YY] = ndgrid(UX,UY);
laser_amplitude = exp(-((XX.^2+YY.^2)/Setup.laserradius^2));
%URX, URY is the real space axis limit, note the non square pixel size
%Initialize with desired amplitude + random phase
mask = double(sqrt(XX.^2+YY.^2)<0.0001);

%Normalize target intensity
for j = 1:LP
    target_amplitude(:,:,j) = target_amplitude(:,:,j)/sqrt(sum(sum(abs(target_amplitude(:,:,j)).^2)));
    %Iteration
end

FieldA = laser_amplitude.*DMDPatternref;

%go back to image plane
[FieldB,psx,psy] = function_lens(FieldA,Setup.ps,Setup.ps,-Setup.f,Setup.lambda);
%Save a copy of the rendering

%Normalize FieldB
FieldB = FieldB/sqrt((sum(abs(FieldB(:).^2))));
EnergyinSignal = sum(sum(abs((1-mask).*FieldB).^2));

fieldNB=FieldB-FieldB;
for j = 1:LP
    %Go from zero to z propoagate the signal only
    fieldZ = function_propagate((1-mask).*FieldB,Setup.lambda,z(j),psy, psx);
    
    %Update the amplitude
    fieldZ = sqrt(EnergyinSignal)*target_amplitude(:,:,j).*exp(1i*angle(fieldZ));
    
    %Go from z to zero propoagate the signal only
    fieldNB =  fieldNB+ function_propagate(fieldZ,Setup.lambda,-z(j),psy, psx);
end
fieldNB = fieldNB/sqrt((sum(abs(fieldNB(:).^2))));

%Add the DC term back at z=0
FieldB = (mask).*FieldB + (1-mask).*fieldNB;
% start in z=0 Go to DMD plane
[FieldA,~,~] = function_lens(FieldB,psx,psy,Setup.f,Setup.lambda);

%Threshold at median value to fill about half the DMD
Amplitude = abs(FieldA);
DMDPattern = double(Amplitude>(median(Amplitude)));

end


function [DMDPattern] = function_GS_Z(Setup, N,z,target_amplitude)
%This an implementation of the Gerchberg–Saxton algorithm 
%to optimize coherent wave modulation patterns.
%Inputs:
%target_amplitude = target 3D intensity distribution
%N = number of iterations
%z = location of target along the z plane
%Outputs:
%DMDPattern = an optimized DMD pattern

UX = Setup.ps*(1:Setup.DMDX); UX = UX-mean(UX);
UY = Setup.ps*(1:Setup.DMDY); UY = UY-mean(UY);
[XX,YY] = ndgrid(UX,UY);
laser_amplitude = exp(-((XX.^2+YY.^2)/Setup.laserradius^2));
FieldB = (laser_amplitude).*exp(1i*2*pi*rand(Setup.DMDX,Setup.DMDY));
[FieldA,psx,psy] = function_lens(FieldB,Setup.ps,Setup.ps,Setup.f,Setup.lambda);
URX = psx*(1:Setup.DMDX); URX = URX-mean(URX);
URY = psy*(1:Setup.DMDY); URY = URY-mean(URY);
[RXX,RYY] = ndgrid(URX,URY);
%URX, URY is the real space axis limit, note the non square pixel size
%Initialize with Desired amplitude and a random phase
FieldZ = (target_amplitude).*exp(1i*2*pi*rand(Setup.DMDX, Setup.DMDY));
mask = double(sqrt(XX.^2+YY.^2)<0.00001);

%Normalize target intensity
target_amplitude = target_amplitude/sqrt((sum(abs(target_amplitude(:).^2))));
%Iteration
for i = 1:N
    
    % start in z=0 Go to DMD plane
    [FieldA,~,~] = function_lens(FieldB,psx,psy,Setup.f,Setup.lambda);
    
    %Threshold at median value to fill about half the DMD
    Amplitude = abs(FieldA);
    DMDPattern = double(Amplitude>(median(Amplitude)));
    FieldA = laser_amplitude.*DMDPattern;
    
    %go back to image plane
    [FieldB,psx,psy] = function_lens(FieldA,Setup.ps,Setup.ps,-Setup.f,Setup.lambda);
    %Save a copy of the rendering
 
    %Normalize FieldB
    FieldB = FieldB/sqrt((sum(abs(FieldB(:).^2))));
    EnergyinSignal = sum(sum(abs((1-mask).*FieldB).^2));
  
    
    %Go from zero to z propoagate the signal only
    fieldZ = function_propagate((1-mask).*FieldB,Setup.lambda,z,psy, psx);
   
    %Update the amplitude
    fieldZ = sqrt(EnergyinSignal)*target_amplitude.*exp(1i*angle(fieldZ));
    
    %Go from z to zero propoagate the signal only
    fieldNB = function_propagate(fieldZ,Setup.lambda,-z,psy, psx);
    
    %Add the DC term back at z=0
    FieldB = (mask).*FieldB + (1-mask).*fieldNB;
    
    %Optional Displays
   
    
end





end


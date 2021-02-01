clear all; close all; clc;

%Setup Properties
Setup.ps = 7.56e-6;         %Pixelsize at DMD
Setup.DMDX = 1080;          % number of pixels on the DMD along one axis
Setup.DMDY = 1920;          % number of pixels on the DMD along the other acis
Setup.laserradius = 200e-4; %Radius of the gaussian laser spot on DMD
Setup.f = 0.2;              % focal length, f, of the telescope lens at an f-f distance between the DMD and the image plane. 
Setup.lambda = 650e-9;      % Wavelength of the laser light source.
N = 10;                     %Number of iterations for the optimization
P = 15;                     %number of coherent patterns being time averaged 
z=[0.05 0.1];               % Depth of the sampling planes, in an array

Zlevels = numel(z);
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

targetIntensity = imread('Target_Images/SmileCat.jpg');
target_amplitude = sqrt(double(mean(targetIntensity,3)));
target_amplitude = target_amplitude/max(target_amplitude(:));
target_amplitude = imresize(target_amplitude', [Setup.DMDX Setup.DMDY]);
target_amplitudeXL = [target_amplitude rot90(target_amplitude,2)];
ta(:,:,1) = imresize(target_amplitudeXL,[Setup.DMDX Setup.DMDY]);
targetIntensity = imread('Target_Images/Cat.jpg');
target_amplitude = sqrt(double(mean(targetIntensity,3)));
target_amplitude = target_amplitude/max(target_amplitude(:));
target_amplitude = imresize(target_amplitude', [Setup.DMDX Setup.DMDY]);
target_amplitudeXL = [target_amplitude rot90(target_amplitude,2)];
ta(:,:,2) = imresize(target_amplitudeXL,[Setup.DMDX Setup.DMDY]);

for kk = 1:Zlevels
    ta(:,:,kk) = squeeze(ta(:,:,kk))/sqrt((sum(sum(sum(abs(ta(:,:,kk).^2))))));
end

for k = 1:P
DMDPatterns(:,:,k) = function_GS_Z(Setup, 1,z(1),ta(:,:,1));
end

[VolumeImages] = function_Rendering(Setup, z,DMDPatterns);

for i = 1:N
VolumeImageAV = mean(VolumeImages,4);

%Display
f = figure(1);
subplot(2,3,1)
imagesc(URX,URY,ta(:,:,1)); colormap gray; axis off;
title('Target Amplitude 1')
subplot(2,3,2)
imagesc(URX,URY,ta(:,:,2)); colormap gray; axis off;
title('Target Amplitude 2')
subplot(2,3,4)
imagesc(URX,URY,VolumeImageAV(:,:,1)); colormap gray; axis off;
title('DCGH Rendering 1')
subplot(2,3,5)
imagesc(URX,URY,VolumeImageAV(:,:,2)); colormap gray; axis off;
title('DCGH Rendering 2')
subplot(2,3,3)
scatter(i,function_accuracy(VolumeImageAV,ta.^2)); hold on;
title('Accuracy per iteration')
ylabel('Accuracy [A.U.]')
xlabel('Iteration [A.U.]')

%New target amplitude computation 
for k = 1:P
newa = ta-ta;
for kk = 1:Zlevels
   baseintensity = ta(:,:,kk).^2;
for q = 1:P
    if q ~= k
    baseintensity=baseintensity-VolumeImages(:,:,kk,q)/P;
    end
end
newa(:,:,kk) = sqrt(max(baseintensity,0));
end

%Use target amplitude to update DMD patterns
DMDPatterns(:,:,k) = function_GS_Z3D_inc(Setup, DMDPatterns(:,:,k),z,newa);

end

%Render images from DMD patterns
[VolumeImages] = function_Rendering(Setup, z,DMDPatterns);


drawnow
end

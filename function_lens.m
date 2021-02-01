function [field2,psX,psY] = function_lens(field1,psx,psy,f,lambda)
%This function simulates the propagation of light through an f-f system : 
% a portion of free space of thickness f, a lens of focal length f>0, and a
% second portion of free space of thickness f.
%Inputs:
%field1 = the input complex field
%psy and psy = pixel size
%f = focal distance of the lens
%lambda = wavelength
%Outputs:
%field2 = the output field in Fourier space 
%psX and psY = new pixel size  

if f<0
field2=ifftshift(ifft2(fftshift(field1)));
else
field2=ifftshift(fft2(fftshift(field1)));
end
[NX,NY] = size(field1);
psX = abs(f*lambda/(NX*psx));
psY = abs(f*lambda/(NY*psy));
end
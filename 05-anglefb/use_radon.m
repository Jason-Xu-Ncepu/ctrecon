clear; clc

N = 256;
focal = 250;
delta_detector = 0.25;

I = phantom(N);
P = fanbeam(I, focal, 'FanSensorGeometry', 'arc', 'FanSensorSpacing', delta_detector);

figure; 
subplot(121), imshow(I, []), title('Shepp-Logan');
subplot(122), imagesc(P), colormap(gray), colorbar, title('Fan-beam sinogram');

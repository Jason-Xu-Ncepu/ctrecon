clear; clc

N = 256;
I = phantom(N);

theta = 0 : 179; % 投影角度
P = radon(I, theta);

figure; imshow(I, []), title('256 Shepp Logan');
figure; imagesc(P), colormap('gray'), colorbar, title('180 degrees projections');
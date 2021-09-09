clear; clc

%% 生成投影数据 %%
N = 256;
I = phantom(N);
theta = 0 : 179;
P = radon(I, theta);

%% 反 radon 变换 %%
rec = iradon(P, theta, 'None'); % 无滤波器，直接反投影
figure; imshow(I, []); title('Original Shepp-Logan');
figure; imshow(rec, []); title('Reconstructed');
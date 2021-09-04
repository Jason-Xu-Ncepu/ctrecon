clear; clc

%% 生成投影数据 %%
N = 256;
I = phantom(N);
theta = 0 : 179;
P = radon(I, theta);

%% 滤波反投影 %%
rec_BP = iradon(P, theta, 'None'); % 直接反投影
rec_RL = iradon(P, theta); % R-L 滤波器
rec_SL = iradon(P, theta, 'linear', 'Shepp-Logan'); % S-L 滤波器

figure; subplot(221); imshow(I, []); title('Origin');
subplot(222); imshow(rec_BP, []); title('BP');
subplot(223); imshow(rec_RL, []); title('RL');
subplot(224); imshow(rec_SL, []); title('SL');
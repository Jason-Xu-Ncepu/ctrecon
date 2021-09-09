clear; clc

N = 256;
N_d = 256;
I = phantom(N);
delta = pi / 180;
theta = 0 : 179;
len_theta = length(theta);

P = radon(I, theta);
[mm, nn] = size(P);
e = floor((mm - N - 1) / 2 + 1) + 1;
P = P(e : N + e - 1, :);
P1 = reshape(P, N, len_theta);

recon = dbpRecon(len_theta, N, P1, delta);

figure; imshow(I, []), title('256 Shepp Logan Original');
figure; imshow(recon, []), title('256 Shepp Logan DBP');


function [recon] = dbpRecon(len_theta, N, R1, delta)
  % len_theta: 投影角度数量；N：图像大小和探测器数量；R1：投影数据（N*len_theta）；delta：角度增量（RAD）
  recon = zeros(N);
  
  for m = 1 : len_theta
    pm = R1(:, m);
    Cm = (N / 2) * (1 - cos((m - 1) * delta) - sin((m - 1) * delta));

    for k1 = 1 : N
      for k2 = 1 : N
        Xrm = Cm + (k2 - 1) * cos((m - 1) * delta) + (k1 - 1) * sin((m - 1) * delta);
        n = floor(Xrm); % 射束编号
        t = Xrm - n; % 小数部分，需要插值
        n = max(1, n); n = min(n, N - 1);
        p = (1 - t) * pm(n) + t * pm(n + 1);
        recon(N + 1 - k1, k2) = recon(N + 1 - k1, k2) + p; % 填充反投影结果
      end
    end
  end
end
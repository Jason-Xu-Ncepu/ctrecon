clear; clc

N = 256;
I = phantom(N);
delta = pi / 180;
theta = 0 : 179;
len_theta = length(theta);
d = 1;

figure;

P = radon(I, theta);
subplot(121), imagesc(P), colormap(gray), title('I sinogram');

[mm, nn] = size(P);
e = floor((mm - N - 1) / 2 + 1) + 1;
P = P(e : N + e - 1, :); % 截取中心N点投影数据
subplot(122), imagesc(P), colormap(gray), title('Clipped I sinogram');

% fakeN = 368;
recon_RL = recon(len_theta, N, P, delta, RLFilter(N, d));
recon_SL = recon(len_theta, N, P, delta, SLFilter(N, d));

figure;
subplot(131), imshow(I, []), title('256 Shepp-Logan');
subplot(132), imshow(recon_RL, []), title('recon RL');
subplot(133), imshow(recon_SL, []), title('recon SL');

function [res] = RLFilter(N, d)
  res = zeros(1, N);
  for k1 = 1 : N
    if mod(k1 - N / 2 - 1, 2) ~= 0
      res(k1) = -1 / (pi * pi * ((k1 - N / 2 - 1) * d) ^ 2);
    end
  end
  res(N / 2 + 1) = 1 / (4 * d ^ 2);
end

function [res] = SLFilter(N, d)
  res = zeros(1, N);
  for k1 = 1 : N
    res(k1) = -2 / (pi ^ 2 * d ^ 2 * (4 * (k1 - N / 2 - 1) ^ 2 - 1));
  end
end

function [recon_] = recon(len_theta, N, R1, delta, filter)
  recon_ = zeros(N);
  for m = 1 : len_theta
    pm = R1(:, m);
    pm_RL = conv(filter, pm, 'same'); % 滤波
    Cm = (N / 2) * (1 - cos((m - 1) * delta) - sin((m - 1) * delta));
    for k1 = 1 : N
      for k2 = 1 : N
        Xrm = Cm + (k2 - 1) * cos((m - 1) * delta) + (k1 - 1) * sin((m - 1) * delta);
        n = floor(Xrm); % 射束编号
        t = Xrm - n; % 小数部分，需要插值
        n = max(1, n); n = min(n, N - 1);
        p = (1 - t) * pm_RL(n) + t * pm_RL(n + 1);
        recon_(N + 1 - k1, k2) = recon_(N + 1 - k1, k2) + p; % 填充反投影结果
      end
    end
  end
end

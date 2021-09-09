clear; clc

N = 256;
theta = 0 : 179;
I = phantom(N);
N_d = 2 * ceil(norm(size(I) - floor((size(I) - 1) / 2) - 1)) + 3; % 探测器数量

P = parallelBeamProjection(theta, N, N_d);

figure;
imshow(I, []), title('Shepp-Logan');
figure;
imagesc(P), colormap(gray), colorbar, title('180 degrees projections');

function [P] = parallelBeamProjection(theta, N, N_d)
  % theta: 投影角度（DEG）；N：图像大小；N_d：探测器数量
  % P: 投影数据（N_d * len(theta))
  shep = [
    % \rho a b x0 y0 \phi
    1 0.69  0.92  0 0 0
    -0.8  0.6624  0.8740  0 -0.0184 0
    -0.2  0.11  0.31  0.22  0 -18
    -0.2  0.16  0.41  -0.22 0 18
    0.1 0.21  0.25  0 0.35  0
    0.1 0.046 0.046 0 0.1 0
    0.1 0.046 0.046 0 -0.1  0
    0.1 0.046 0.023 -0.08 0.605 0
    0.1 0.023 0.023 0 -0.606  0
    0.1 0.023 0.046 0.06  -0.605  0
  ]; % Shepp-Logan 参数
  theta = theta * pi / 180;
  len_theta = length(theta);
  P = zeros(N_d, len_theta); % 投影数据结果

  rho = shep(:, 1).';
  ae = 0.5 * N * shep(:, 2).';
  be = 0.5 * N * shep(:, 3).';
  xe = 0.5 * N * shep(:, 4).';
  ye = 0.5 * N * shep(:, 5).';
  alpha = shep(:, 6).' * pi / 180;

  detectors = -(N_d - 1) / 2 : (N_d - 1) / 2;

  for k1 = 1 :len_theta
    P_theta = zeros(1, N_d);
    for k2 = 1 : max(size(xe))
      a = (ae(k2) * cos(theta(k1) - alpha(k2))) ^ 2 + (be(k2) * sin(theta(k1) - alpha(k2))) .^ 2;
      temp = a - (detectors - xe(k2) * cos(theta(k1)) - ye(k2) * sin(theta(k1))) .^ 2;
      ind = temp > 0;
      P_theta(ind) = P_theta(ind) + rho(k2) * (2 * ae(k2) * be(k2) * sqrt(temp(ind))) ./ a;
    end
      P(:, k1) = P_theta.';
  end
end
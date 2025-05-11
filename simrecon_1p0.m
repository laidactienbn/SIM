clear all; clc;
%% read all pngs and do FT

imraw = im2double(imread('unprocessed/imraw.png'));
imconv = im2double(imread('unprocessed/imconv.png'));

im_vert_0pi_nl = im2double(imread('unprocessed/im_vert_0pi_nl.png'));
im_vert_0p4pi_nl = im2double(imread('unprocessed/im_vert_0p4pi_nl.png'));
im_vert_0p8pi_nl = im2double(imread('unprocessed/im_vert_0p8pi_nl.png'));
im_vert_1p2pi_nl = im2double(imread('unprocessed/im_vert_1p2pi_nl.png'));
im_vert_1p6pi_nl = im2double(imread('unprocessed/im_vert_1p6pi_nl.png'));

im_horz_0pi_nl = im2double(imread('unprocessed/im_horz_0pi_nl.png'));
im_horz_0p4pi_nl = im2double(imread('unprocessed/im_horz_0p4pi_nl.png'));
im_horz_0p8pi_nl = im2double(imread('unprocessed/im_horz_0p8pi_nl.png'));
im_horz_1p2pi_nl = im2double(imread('unprocessed/im_horz_1p2pi_nl.png'));
im_horz_1p6pi_nl = im2double(imread('unprocessed/im_horz_1p6pi_nl.png'));

im_diag_l45_0pi_nl = im2double(imread('unprocessed/im_diag_l45_0pi_nl.png'));
im_diag_l45_0p4pi_nl = im2double(imread('unprocessed/im_diag_l45_0p4pi_nl.png'));
im_diag_l45_0p8pi_nl = im2double(imread('unprocessed/im_diag_l45_0p8pi_nl.png'));
im_diag_l45_1p2pi_nl = im2double(imread('unprocessed/im_diag_l45_1p2pi_nl.png'));
im_diag_l45_1p6pi_nl = im2double(imread('unprocessed/im_diag_l45_1p6pi_nl.png'));

im_diag_r45_0pi_nl = im2double(imread('unprocessed/im_diag_r45_0pi_nl.png'));
im_diag_r45_0p4pi_nl = im2double(imread('unprocessed/im_diag_r45_0p4pi_nl.png'));
im_diag_r45_0p8pi_nl = im2double(imread('unprocessed/im_diag_r45_0p8pi_nl.png'));
im_diag_r45_1p2pi_nl = im2double(imread('unprocessed/im_diag_r45_1p2pi_nl.png'));
im_diag_r45_1p6pi_nl = im2double(imread('unprocessed/im_diag_r45_1p6pi_nl.png'));

psfraw = im2double(imread('unprocessed/psfraw.png'));

psfraw_vert_0pi_nl = im2double(imread('unprocessed/psfraw_vert_0pi_nl.png'));
psfraw_vert_0p4pi_nl = im2double(imread('unprocessed/psfraw_vert_0p4pi_nl.png'));
psfraw_vert_0p8pi_nl = im2double(imread('unprocessed/psfraw_vert_0p8pi_nl.png'));
psfraw_vert_1p2pi_nl = im2double(imread('unprocessed/psfraw_vert_1p2pi_nl.png'));
psfraw_vert_1p6pi_nl = im2double(imread('unprocessed/psfraw_vert_1p6pi_nl.png'));

psfraw_horz_0pi_nl = im2double(imread('unprocessed/psfraw_horz_0pi_nl.png'));
psfraw_horz_0p4pi_nl = im2double(imread('unprocessed/psfraw_horz_0p4pi_nl.png'));
psfraw_horz_0p8pi_nl = im2double(imread('unprocessed/psfraw_horz_0p8pi_nl.png'));
psfraw_horz_1p2pi_nl = im2double(imread('unprocessed/psfraw_horz_1p2pi_nl.png'));
psfraw_horz_1p6pi_nl = im2double(imread('unprocessed/psfraw_horz_1p6pi_nl.png'));

psfraw_diag_l45_0pi_nl = im2double(imread('unprocessed/psfraw_diag_l45_0pi_nl.png'));
psfraw_diag_l45_0p4pi_nl = im2double(imread('unprocessed/psfraw_diag_l45_0p4pi_nl.png'));
psfraw_diag_l45_0p8pi_nl = im2double(imread('unprocessed/psfraw_diag_l45_0p8pi_nl.png'));
psfraw_diag_l45_1p2pi_nl = im2double(imread('unprocessed/psfraw_diag_l45_1p2pi_nl.png'));
psfraw_diag_l45_1p6pi_nl = im2double(imread('unprocessed/psfraw_diag_l45_1p6pi_nl.png'));

psfraw_diag_r45_0pi_nl = im2double(imread('unprocessed/psfraw_diag_r45_0pi_nl.png'));
psfraw_diag_r45_0p4pi_nl = im2double(imread('unprocessed/psfraw_diag_r45_0p4pi_nl.png'));
psfraw_diag_r45_0p8pi_nl = im2double(imread('unprocessed/psfraw_diag_r45_0p8pi_nl.png'));
psfraw_diag_r45_1p2pi_nl = im2double(imread('unprocessed/psfraw_diag_r45_1p2pi_nl.png'));
psfraw_diag_r45_1p6pi_nl = im2double(imread('unprocessed/psfraw_diag_r45_1p6pi_nl.png'));

mat_vert_0pi = im2double(imread('unprocessed/mat_vert_0pi.png'));
mat_horz_0pi = im2double(imread('unprocessed/mat_horz_0pi.png'));
mat_diag_l45_0pi = im2double(imread('unprocessed/mat_diag_l45_0pi.png'));
mat_diag_r45_0pi = im2double(imread('unprocessed/mat_diag_r45_0pi.png'));



sz = size(im_vert_0pi_nl, 1);

ft_im_vert_0pi_nl = fftshift(fft2(im_vert_0pi_nl));
ft_im_vert_0p4pi_nl = fftshift(fft2(im_vert_0p4pi_nl));
ft_im_vert_0p8pi_nl = fftshift(fft2(im_vert_0p8pi_nl));
ft_im_vert_1p2pi_nl = fftshift(fft2(im_vert_1p2pi_nl));
ft_im_vert_1p6pi_nl = fftshift(fft2(im_vert_1p6pi_nl));

ft_im_horz_0pi_nl = fftshift(fft2(im_horz_0pi_nl));
ft_im_horz_0p4pi_nl = fftshift(fft2(im_horz_0p4pi_nl));
ft_im_horz_0p8pi_nl = fftshift(fft2(im_horz_0p8pi_nl));
ft_im_horz_1p2pi_nl = fftshift(fft2(im_horz_1p2pi_nl));
ft_im_horz_1p6pi_nl = fftshift(fft2(im_horz_1p6pi_nl));

ft_im_diag_l45_0pi_nl = fftshift(fft2(im_diag_l45_0pi_nl));
ft_im_diag_l45_0p4pi_nl = fftshift(fft2(im_diag_l45_0p4pi_nl));
ft_im_diag_l45_0p8pi_nl = fftshift(fft2(im_diag_l45_0p8pi_nl));
ft_im_diag_l45_1p2pi_nl = fftshift(fft2(im_diag_l45_1p2pi_nl));
ft_im_diag_l45_1p6pi_nl = fftshift(fft2(im_diag_l45_1p6pi_nl));

ft_im_diag_r45_0pi_nl = fftshift(fft2(im_diag_r45_0pi_nl));
ft_im_diag_r45_0p4pi_nl = fftshift(fft2(im_diag_r45_0p4pi_nl));
ft_im_diag_r45_0p8pi_nl = fftshift(fft2(im_diag_r45_0p8pi_nl));
ft_im_diag_r45_1p2pi_nl = fftshift(fft2(im_diag_r45_1p2pi_nl));
ft_im_diag_r45_1p6pi_nl = fftshift(fft2(im_diag_r45_1p6pi_nl));

otfraw = fftshift(fft2(psfraw));

ft_psfraw_vert_0pi_nl = fftshift(fft2(psfraw_vert_0pi_nl));
ft_psfraw_vert_0p4pi_nl = fftshift(fft2(psfraw_vert_0p4pi_nl));
ft_psfraw_vert_0p8pi_nl = fftshift(fft2(psfraw_vert_0p8pi_nl));
ft_psfraw_vert_1p2pi_nl = fftshift(fft2(psfraw_vert_1p2pi_nl));
ft_psfraw_vert_1p6pi_nl = fftshift(fft2(psfraw_vert_1p6pi_nl));

ft_psfraw_horz_0pi_nl = fftshift(fft2(psfraw_horz_0pi_nl));
ft_psfraw_horz_0p4pi_nl = fftshift(fft2(psfraw_horz_0p4pi_nl));
ft_psfraw_horz_0p8pi_nl = fftshift(fft2(psfraw_horz_0p8pi_nl));
ft_psfraw_horz_1p2pi_nl = fftshift(fft2(psfraw_horz_1p2pi_nl));
ft_psfraw_horz_1p6pi_nl = fftshift(fft2(psfraw_horz_1p6pi_nl));

ft_psfraw_diag_l45_0pi_nl = fftshift(fft2(psfraw_diag_l45_0pi_nl));
ft_psfraw_diag_l45_0p4pi_nl = fftshift(fft2(psfraw_diag_l45_0p4pi_nl));
ft_psfraw_diag_l45_0p8pi_nl = fftshift(fft2(psfraw_diag_l45_0p8pi_nl));
ft_psfraw_diag_l45_1p2pi_nl = fftshift(fft2(psfraw_diag_l45_1p2pi_nl));
ft_psfraw_diag_l45_1p6pi_nl = fftshift(fft2(psfraw_diag_l45_1p6pi_nl));

ft_psfraw_diag_r45_0pi_nl = fftshift(fft2(psfraw_diag_r45_0pi_nl));
ft_psfraw_diag_r45_0p4pi_nl = fftshift(fft2(psfraw_diag_r45_0p4pi_nl));
ft_psfraw_diag_r45_0p8pi_nl = fftshift(fft2(psfraw_diag_r45_0p8pi_nl));
ft_psfraw_diag_r45_1p2pi_nl = fftshift(fft2(psfraw_diag_r45_1p2pi_nl));
ft_psfraw_diag_r45_1p6pi_nl = fftshift(fft2(psfraw_diag_r45_1p6pi_nl));

ft_mat_vert_0pi = fftshift(fft2(mat_vert_0pi));
ft_mat_horz_0pi = fftshift(fft2(mat_horz_0pi));
ft_mat_diag_l45_0pi = fftshift(fft2(mat_diag_l45_0pi));
ft_mat_diag_r45_0pi = fftshift(fft2(mat_diag_r45_0pi));


%% reconstruction based on A. Lal, C. Shan and P. Xi, "Structured Illumination Microscopy Image Reconstruction
% Algorithm," in IEEE Journal of Selected Topics in Quantum Electronics, vol. 22, no. 4, pp.
% 50-63, July-Aug. 2016, Art no. 6803414, doi: 10.1109/JSTQE.2016.2521542.

% step 1: solve matrix M
m0 = 1;
m1 = 1;
m2 = 20; % weight factor, set to 1 for now
% m0 = 1;
% m1 = 1;
% m2 = 60; % weight factor, set to 1 for now
phi1 = 0;
phi2 = 2 * pi / 5;
phi3 = 4 * pi / 5;
phi4 = 6 * pi / 5;
phi5 = 8 * pi / 5;

M = [1, - m1 * exp(-1i * phi1) / 2, - m1 * exp(1i * phi1) / 2, - m2 * exp(-2 * 1i * phi1) / 2, - m2 * exp(2 * 1i * phi1) / 2;
    1, - m1 * exp(-1i * phi2) / 2, - m1 * exp(1i * phi2) / 2, - m2 * exp(-2 * 1i * phi2) / 2, - m2 * exp(2 * 1i * phi2) / 2;
    1, - m1 * exp(-1i * phi3) / 2, - m1 * exp(1i * phi3) / 2, - m2 * exp(-2 * 1i * phi3) / 2, - m2 * exp(2 * 1i * phi3) / 2;
    1, - m1 * exp(-1i * phi4) / 2, - m1 * exp(1i * phi4) / 2, - m2 * exp(-2 * 1i * phi4) / 2, - m2 * exp(2 * 1i * phi4) / 2;
    1, - m1 * exp(-1i * phi5) / 2, - m1 * exp(1i * phi5) / 2, - m2 * exp(-2 * 1i * phi5) / 2, - m2 * exp(2 * 1i * phi5) / 2;];
% matrix M, eq(4) in paper
M_inv = inv(M);

D_tilde_vert = {ft_im_vert_0pi_nl;
    ft_im_vert_0p4pi_nl;
    ft_im_vert_0p8pi_nl;
    ft_im_vert_1p2pi_nl;
    ft_im_vert_1p6pi_nl;};

D_tilde_horz = {ft_im_horz_0pi_nl;
    ft_im_horz_0p4pi_nl;
    ft_im_horz_0p8pi_nl;
    ft_im_horz_1p2pi_nl;
    ft_im_horz_1p6pi_nl;};

D_tilde_diag_l45 = {ft_im_diag_l45_0pi_nl;
    ft_im_diag_l45_0p4pi_nl;
    ft_im_diag_l45_0p8pi_nl;
    ft_im_diag_l45_1p2pi_nl;
    ft_im_diag_l45_1p6pi_nl;};

D_tilde_diag_r45 = {ft_im_diag_r45_0pi_nl;
    ft_im_diag_r45_0p4pi_nl;
    ft_im_diag_r45_0p8pi_nl;
    ft_im_diag_r45_1p2pi_nl;
    ft_im_diag_r45_1p6pi_nl;};

S_tilde_H_tilde_vert = {0; 0; 0; 0; 0;};
for i = 1:5;
    for j = 1:5;
        S_tilde_H_tilde_vert{i, 1} =  S_tilde_H_tilde_vert{i, 1} + M_inv(i, j) .* D_tilde_vert{j, 1};
    end;
end;
% equation (5) in paper

Su_tilde_vert = S_tilde_H_tilde_vert;
% assume no noise so no weiner filtering

S_tilde_H_tilde_horz = {0; 0; 0; 0; 0;};
for i = 1:5;
    for j = 1:5;
        S_tilde_H_tilde_horz{i, 1} =  S_tilde_H_tilde_horz{i, 1} + M_inv(i, j) .* D_tilde_horz{j, 1};
    end;
end;
% equation (5) in paper

Su_tilde_horz = S_tilde_H_tilde_horz;
% assume no noise so no weiner filtering

S_tilde_H_tilde_diag_l45 = {0; 0; 0; 0; 0;};
for i = 1:5;
    for j = 1:5;
        S_tilde_H_tilde_diag_l45{i, 1} =  S_tilde_H_tilde_diag_l45{i, 1} + M_inv(i, j) .* D_tilde_diag_l45{j, 1};
    end;
end;
% equation (5) in paper

Su_tilde_diag_l45 = S_tilde_H_tilde_diag_l45;
% assume no noise so no weiner filtering

S_tilde_H_tilde_diag_r45 = {0; 0; 0; 0; 0;};
for i = 1:5;
    for j = 1:5;
        S_tilde_H_tilde_diag_r45{i, 1} =  S_tilde_H_tilde_diag_r45{i, 1} + M_inv(i, j) .* D_tilde_diag_r45{j, 1};
    end;
end;
% equation (5) in paper

Su_tilde_diag_r45 = S_tilde_H_tilde_diag_r45;
% assume no noise so no weiner filtering

% step 2: shift frequency component to true position

%find locations of harmonics in freq space


for i = 1 : 5;
    [val, loc] = max(S_tilde_H_tilde_vert{i, 1}, [], 'all');
    [row, col] = ind2sub(sz, loc);
    harmonic_loc_vert{i, 1} = [row, col];
    [val, loc] = max(S_tilde_H_tilde_horz{i, 1}, [], 'all');
    [row, col] = ind2sub(sz, loc);
    harmonic_loc_horz{i, 1} = [row, col];
    [val, loc] = max(S_tilde_H_tilde_diag_l45{i, 1}, [], 'all');
    [row, col] = ind2sub(sz, loc);
    harmonic_loc_diag_l45{i, 1} = [row, col];
    [val, loc] = max(S_tilde_H_tilde_diag_r45{i, 1}, [], 'all');
    [row, col] = ind2sub(sz, loc);
    harmonic_loc_diag_r45{i, 1} = [row, col];
end

%dist = 70; %filtering distance
[pks, locs, dist, p] = findpeaks(abs(otfraw(sz / 2, :))); % find filtering distance using OTF
cutoff_ratio = 0.7;
dist = dist * cutoff_ratio; % ratio to limit filtering distance, affects OTF shape

for i = 1 : 5;
    for j = 1 : sz;
        for k = 1 : sz;
            if sqrt( (j - sz / 2) ^ 2 + (k - sz / 2) ^ 2) > dist;
                Su_tilde_vert{i, 1}(j, k) = 0;
            end;
            if sqrt( (j - sz / 2) ^ 2 + (k - sz / 2) ^ 2) > dist;
                Su_tilde_horz{i, 1}(j, k) = 0;
            end;
            if sqrt( (j - sz / 2) ^ 2 + (k - sz / 2) ^ 2) > dist;
                Su_tilde_diag_l45{i, 1}(j, k) = 0;
            end;
            if sqrt( (j - sz / 2) ^ 2 + (k - sz / 2) ^ 2) > dist;
                Su_tilde_diag_r45{i, 1}(j, k) = 0;
            end;
        end;
    end;
end;
    % filter out all freqs outside filtering distance

S_shifted_tilde_vert = {Su_tilde_vert{1, 1}; 0; 0; 0; 0;};
S_shifted_tilde_horz = {Su_tilde_horz{1, 1}; 0; 0; 0; 0;};
S_shifted_tilde_diag_l45 = {Su_tilde_diag_l45{1, 1}; 0; 0; 0; 0;};
S_shifted_tilde_diag_r45 = {Su_tilde_diag_r45{1, 1}; 0; 0; 0; 0;};
for j = 2 : 5;
    S_shifted_tilde_vert{j, 1} = fft2(ifft2(Su_tilde_vert{j, 1}) * exp(((-1) ^ (j - 1)) * 1i * 2 * pi * pdist2(harmonic_loc_vert{1, 1}, harmonic_loc_vert{j, 1})));
    S_shifted_tilde_horz{j, 1} = fft2(ifft2(Su_tilde_horz{j, 1}) * exp(((-1) ^ (j - 1)) * 1i * 2 * pi * pdist2(harmonic_loc_horz{1, 1}, harmonic_loc_horz{j, 1})));
    S_shifted_tilde_diag_r45{j, 1} = fft2(ifft2(Su_tilde_diag_l45{j, 1}) * exp(((-1) ^ (j - 1)) * 1i * 2 * pi * pdist2(harmonic_loc_diag_l45{1, 1}, harmonic_loc_diag_l45{j, 1})));
    S_shifted_tilde_diag_l45{j, 1} = fft2(ifft2(Su_tilde_diag_r45{j, 1}) * exp(((-1) ^ (j - 1)) * 1i * 2 * pi * pdist2(harmonic_loc_diag_r45{1, 1}, harmonic_loc_diag_r45{j, 1})));
end;
%equation 14, 15 in paper

%% step 3: shift freq components so they line up and are all centered at 0
for j = 2 : 5;
    S_shifted_tilde_moved_vert{j, 1} = circshift(S_shifted_tilde_vert{j , 1}, harmonic_loc_vert{1, 1} - harmonic_loc_vert{j, 1});
    S_shifted_tilde_moved_horz{j, 1} = circshift(S_shifted_tilde_horz{j , 1}, harmonic_loc_horz{1, 1} - harmonic_loc_horz{j, 1});
    S_shifted_tilde_moved_diag_l45{j, 1} = circshift(S_shifted_tilde_diag_l45{j , 1}, (harmonic_loc_diag_l45{1, 1} - harmonic_loc_diag_l45{j, 1}) .* [1, -1]);
    S_shifted_tilde_moved_diag_r45{j, 1} = circshift(S_shifted_tilde_diag_r45{j , 1}, (harmonic_loc_diag_r45{1, 1} - harmonic_loc_diag_r45{j, 1}) .* [1, -1]);
end;

%% step 4: fuse freq components!


ft_sim = m0 * (((S_shifted_tilde_vert{1, 1} + S_shifted_tilde_horz{1, 1} + S_shifted_tilde_diag_l45{1, 1} + S_shifted_tilde_diag_r45{1, 1})));


for i = 2:3 ;
    for j = 1:sz;
        for k = 1:sz
            ft_sim(j, k) = ft_sim(j, k) + m1^2 * S_shifted_tilde_moved_vert{i, 1}(j, k) + m1^2 * S_shifted_tilde_moved_horz{i, 1}(j, k) + m1^2 * S_shifted_tilde_moved_diag_l45{i, 1}(j, k) + m1^2 * S_shifted_tilde_moved_diag_r45{i, 1}(j, k);
        end
    end
end
for i = 4:5 ;
    for j = 1:sz;
        for k = 1:sz
            ft_sim(j, k) = ft_sim(j, k) + m2^2 * S_shifted_tilde_moved_vert{i, 1}(j, k) + m2^2 * S_shifted_tilde_moved_horz{i, 1}(j, k) + m2^2 * S_shifted_tilde_moved_diag_l45{i, 1}(j, k) + m2^2 * S_shifted_tilde_moved_diag_r45{i, 1}(j, k);
            
        end
    end
end

sim = ifft2(ft_sim);
ft_im_conv = fftshift(fft2(imconv));
ft_raw = fftshift(fft2(imraw));

%Freq
figure;
% Show first image on the left
subplot(1, 3, 1);
imshow(0.1*log(ft_raw));
title('Raw');

% Show first image on the left
subplot(1, 3, 2);
imshow(0.1*log(ft_im_conv));
title('Direct');

% Show second image on the right
subplot(1, 3, 3);
imshow(0.1*log(ft_sim));
title('SIM');

%Real
figure;
% Show first image on the left
subplot(1, 3, 1);
imshow(imraw);
title('Raw');

% Show first image on the left
subplot(1, 3, 2);
imshow(imconv);
title('Direct');

% Show second image on the right
subplot(1, 3, 3);
imshow(sim);
title('SIM');
% imwrite(imraw,'imraw.png');
% imwrite(imconv, 'imconv.png');
% imwrite(sim, 'sim.png');
% imwrite(ft_raw, 'ft_raw.png');
% imwrite(ft_im_conv, 'ft_im_conv.png');
% imwrite(ft_sim, 'ft_sim.png');

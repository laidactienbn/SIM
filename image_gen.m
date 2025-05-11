%% generate images for SIM

clear; clc;

% fringe generation
% since psfraw FWHM of SLAM1 is ~500nm and minimum pixel resolution (step
% size) is 100nm, assume image of 500 pixels x 0.1um step size = 50um FOV,
% assume stripes are 1um or 10 pixels wide

sz = 500; % pixels per rows, cols, must be multiple of 10

% approximate appropriate modulation values for a sine wave
stripepxnum = 10; %number of pixels per stripe
x = linspace(0, 2 * pi, stripepxnum);
modvals = (-cos(x) + 1) / 2;
% modvals = horzcat(modvals,modvals)
% plot(modvals)

%% vertical stripes
stripenum = floor(sz / stripepxnum);
mat = modvals;
for i = 1 : stripenum - 1
    mat = horzcat(mat, modvals);
end;

mat = repmat(mat, sz, 1);
mat_vert_0pi = mat;

% shift stripe phase 4 times by intervals of 2pi/5
shiftpxnum = stripepxnum / 5;

%2pi/5
mat_vert_0p4pi = circshift(mat, shiftpxnum, 2);
%4pi/5
mat_vert_0p8pi = circshift(mat, shiftpxnum * 2, 2);
%6pi/5
mat_vert_1p2pi = circshift(mat, shiftpxnum * 3, 2);
%8pi/5
mat_vert_1p6pi = circshift(mat, shiftpxnum * 4, 2);

%% horizontal stripes
mat_horz_0pi = mat_vert_0pi.';
mat_horz_0p4pi = mat_vert_0p4pi.';
mat_horz_0p8pi = mat_vert_0p8pi.';
mat_horz_1p2pi = mat_vert_1p2pi.';
mat_horz_1p6pi = mat_vert_1p6pi.';

%% diagonal stripes left leaning 45 deg
%0pi
for i = 1:sz;
    mat_diag_l45_0pi(i, :) = circshift(mat_vert_0pi(i, :), i, 2);
end
%2pi/5
for i = 1:sz;
    mat_diag_l45_0p4pi(i, :) = circshift(mat_vert_0p4pi(i, :), i, 2);
end
%4pi/5
for i = 1:sz;
    mat_diag_l45_0p8pi(i, :) = circshift(mat_vert_0p8pi(i, :), i, 2);
end
%6pi/5
for i = 1:sz;
    mat_diag_l45_1p2pi(i, :) = circshift(mat_vert_1p2pi(i, :), i, 2);
end
%8pi/5
for i = 1:sz;
    mat_diag_l45_1p6pi(i, :) = circshift(mat_vert_1p6pi(i, :), i, 2);
end

%% diagonal stripes right leaning 45 deg 
 mat_diag_r45_0pi = rot90(mat_diag_l45_0pi);
 mat_diag_r45_0p4pi = rot90(mat_diag_l45_0p4pi);
 mat_diag_r45_0p8pi = rot90(mat_diag_l45_0p8pi);
 mat_diag_r45_1p2pi = rot90(mat_diag_l45_1p2pi);
 mat_diag_r45_1p6pi = rot90(mat_diag_l45_1p6pi);

 %% due to 2p nonlinear effect (emmission I = excitation I ^ 2), square all components of mat for nonlinear SIM
 mat_vert_0pi_nl = mat_vert_0pi .^ 2;
 mat_vert_0p4pi_nl = mat_vert_0p4pi .^ 2;
 mat_vert_0p8pi_nl = mat_vert_0p8pi .^ 2;
 mat_vert_1p2pi_nl = mat_vert_1p2pi .^ 2;
 mat_vert_1p6pi_nl = mat_vert_1p6pi .^ 2;

 mat_horz_0pi_nl = mat_horz_0pi .^ 2;
 mat_horz_0p4pi_nl = mat_horz_0p4pi .^ 2;
 mat_horz_0p8pi_nl = mat_horz_0p8pi .^ 2;
 mat_horz_1p2pi_nl = mat_horz_1p2pi .^ 2;
 mat_horz_1p6pi_nl = mat_horz_1p6pi .^ 2;

 mat_diag_l45_0pi_nl = mat_diag_l45_0pi .^ 2;
 mat_diag_l45_0p4pi_nl = mat_diag_l45_0p4pi .^ 2;
 mat_diag_l45_0p8pi_nl = mat_diag_l45_0p8pi .^ 2;
 mat_diag_l45_1p2pi_nl = mat_diag_l45_1p2pi .^ 2;
 mat_diag_l45_1p6pi_nl = mat_diag_l45_1p6pi .^ 2;

 mat_diag_r45_0pi_nl = mat_diag_r45_0pi .^ 2;
 mat_diag_r45_0p4pi_nl = mat_diag_r45_0p4pi .^ 2;
 mat_diag_r45_0p8pi_nl = mat_diag_r45_0p8pi .^ 2;
 mat_diag_r45_1p2pi_nl = mat_diag_r45_1p2pi .^ 2;
 mat_diag_r45_1p6pi_nl = mat_diag_r45_1p6pi .^ 2;

 %% TEST
 % n = 1;
 % mat_vert_0pi_nl = n+mat_vert_0pi .^ 2;
 % mat_vert_0p4pi_nl = n+mat_vert_0p4pi .^ 2;
 % mat_vert_0p8pi_nl = n+mat_vert_0p8pi .^ 2;
 % mat_vert_1p2pi_nl = n+mat_vert_1p2pi .^ 2;
 % mat_vert_1p6pi_nl = n+mat_vert_1p6pi .^ 2;
 % 
 % mat_horz_0pi_nl = n+mat_horz_0pi .^ 2;
 % mat_horz_0p4pi_nl = n+mat_horz_0p4pi .^ 2;
 % mat_horz_0p8pi_nl = n+mat_horz_0p8pi .^ 2;
 % mat_horz_1p2pi_nl = n+mat_horz_1p2pi .^ 2;
 % mat_horz_1p6pi_nl = n+mat_horz_1p6pi .^ 2;
 % 
 % mat_diag_l45_0pi_nl = n+mat_diag_l45_0pi .^ 2;
 % mat_diag_l45_0p4pi_nl = n+mat_diag_l45_0p4pi .^ 2;
 % mat_diag_l45_0p8pi_nl = n+mat_diag_l45_0p8pi .^ 2;
 % mat_diag_l45_1p2pi_nl = n+mat_diag_l45_1p2pi .^ 2;
 % mat_diag_l45_1p6pi_nl = n+mat_diag_l45_1p6pi .^ 2;
 % 
 % mat_diag_r45_0pi_nl = n+mat_diag_r45_0pi .^ 2;
 % mat_diag_r45_0p4pi_nl = n+mat_diag_r45_0p4pi .^ 2;
 % mat_diag_r45_0p8pi_nl = n+mat_diag_r45_0p8pi .^ 2;
 % mat_diag_r45_1p2pi_nl = n+mat_diag_r45_1p2pi .^ 2;
 % mat_diag_r45_1p6pi_nl = n+mat_diag_r45_1p6pi .^ 2;

 % mat_vert_0pi_nl = mat_vert_0pi  ;
 % mat_vert_0p4pi_nl = mat_vert_0p4pi  ;
 % mat_vert_0p8pi_nl = mat_vert_0p8pi  ;
 % mat_vert_1p2pi_nl = mat_vert_1p2pi  ;
 % mat_vert_1p6pi_nl = mat_vert_1p6pi  ;
 % 
 % mat_horz_0pi_nl = mat_horz_0pi  ;
 % mat_horz_0p4pi_nl = mat_horz_0p4pi  ;
 % mat_horz_0p8pi_nl = mat_horz_0p8pi  ;
 % mat_horz_1p2pi_nl = mat_horz_1p2pi  ;
 % mat_horz_1p6pi_nl = mat_horz_1p6pi  ;
 % 
 % mat_diag_l45_0pi_nl = mat_diag_l45_0pi  ;
 % mat_diag_l45_0p4pi_nl = mat_diag_l45_0p4pi  ;
 % mat_diag_l45_0p8pi_nl = mat_diag_l45_0p8pi  ;
 % mat_diag_l45_1p2pi_nl = mat_diag_l45_1p2pi  ;
 % mat_diag_l45_1p6pi_nl = mat_diag_l45_1p6pi  ;
 % 
 % mat_diag_r45_0pi_nl = mat_diag_r45_0pi  ;
 % mat_diag_r45_0p4pi_nl = mat_diag_r45_0p4pi  ;
 % mat_diag_r45_0p8pi_nl = mat_diag_r45_0p8pi  ;
 % mat_diag_r45_1p2pi_nl = mat_diag_r45_1p2pi  ;
 % mat_diag_r45_1p6pi_nl = mat_diag_r45_1p6pi  ;


 %% make SIM images

 %load image, must be same size as sz
%  imraw = im2gray(imread('peacock.jpg'));
%  imraw = imraw(230:729, 120:619);
 imraw = im2gray(imread('strands.png'));
 %imraw = im2gray(imread('siemens.png'));
%  imraw = im2gray(imread('unprocessed/psfraw.tiff'));
 %convolve with gaussian kernel
 sigma = 2; % FWHM ~ 2.5 * sigma
 imconv = im2double(imgaussfilt(imraw,sigma));
 %imshow(imconv)

 % multiply by stripes
 im_vert_0pi_nl = imconv .* mat_vert_0pi_nl;
 im_vert_0p4pi_nl = imconv .* mat_vert_0p4pi_nl;
 im_vert_0p8pi_nl = imconv .* mat_vert_0p8pi_nl;
 im_vert_1p2pi_nl = imconv .* mat_vert_1p2pi_nl;
 im_vert_1p6pi_nl = imconv .* mat_vert_1p6pi_nl;

 im_horz_0pi_nl = imconv .* mat_horz_0pi_nl;
 im_horz_0p4pi_nl = imconv .* mat_horz_0p4pi_nl;
 im_horz_0p8pi_nl = imconv .* mat_horz_0p8pi_nl;
 im_horz_1p2pi_nl = imconv .* mat_horz_1p2pi_nl;
 im_horz_1p6pi_nl = imconv .* mat_horz_1p6pi_nl;

 im_diag_l45_0pi_nl = imconv .* mat_diag_l45_0pi_nl;
 im_diag_l45_0p4pi_nl = imconv .* mat_diag_l45_0p4pi_nl;
 im_diag_l45_0p8pi_nl = imconv .* mat_diag_l45_0p8pi_nl;
 im_diag_l45_1p2pi_nl = imconv .* mat_diag_l45_1p2pi_nl;
 im_diag_l45_1p6pi_nl = imconv .* mat_diag_l45_1p6pi_nl;

 im_diag_r45_0pi_nl = imconv .* mat_diag_r45_0pi_nl;
 im_diag_r45_0p4pi_nl = imconv .* mat_diag_r45_0p4pi_nl;
 im_diag_r45_0p8pi_nl = imconv .* mat_diag_r45_0p8pi_nl;
 im_diag_r45_1p2pi_nl = imconv .* mat_diag_r45_1p2pi_nl;
 im_diag_r45_1p6pi_nl = imconv .* mat_diag_r45_1p6pi_nl;

 %% generate psf

 pointsource = zeros(501, 501);
 pointsource(251, 251) = 1;
 psfraw = imgaussfilt(pointsource,sigma);
 psfraw = imresize(psfraw, 500 / 501);

 % multiply by stripes
 psfraw_vert_0pi_nl = psfraw .* mat_vert_0pi_nl;
 psfraw_vert_0p4pi_nl = psfraw .* mat_vert_0p4pi_nl;
 psfraw_vert_0p8pi_nl = psfraw .* mat_vert_0p8pi_nl;
 psfraw_vert_1p2pi_nl = psfraw .* mat_vert_1p2pi_nl;
 psfraw_vert_1p6pi_nl = psfraw .* mat_vert_1p6pi_nl;

 psfraw_horz_0pi_nl = psfraw .* mat_horz_0pi_nl;
 psfraw_horz_0p4pi_nl = psfraw .* mat_horz_0p4pi_nl;
 psfraw_horz_0p8pi_nl = psfraw .* mat_horz_0p8pi_nl;
 psfraw_horz_1p2pi_nl = psfraw .* mat_horz_1p2pi_nl;
 psfraw_horz_1p6pi_nl = psfraw .* mat_horz_1p6pi_nl;

 psfraw_diag_l45_0pi_nl = psfraw .* mat_diag_l45_0pi_nl;
 psfraw_diag_l45_0p4pi_nl = psfraw .* mat_diag_l45_0p4pi_nl;
 psfraw_diag_l45_0p8pi_nl = psfraw .* mat_diag_l45_0p8pi_nl;
 psfraw_diag_l45_1p2pi_nl = psfraw .* mat_diag_l45_1p2pi_nl;
 psfraw_diag_l45_1p6pi_nl = psfraw .* mat_diag_l45_1p6pi_nl;

 psfraw_diag_r45_0pi_nl = psfraw .* mat_diag_r45_0pi_nl;
 psfraw_diag_r45_0p4pi_nl = psfraw .* mat_diag_r45_0p4pi_nl;
 psfraw_diag_r45_0p8pi_nl = psfraw .* mat_diag_r45_0p8pi_nl;
 psfraw_diag_r45_1p2pi_nl = psfraw .* mat_diag_r45_1p2pi_nl;
 psfraw_diag_r45_1p6pi_nl = psfraw .* mat_diag_r45_1p6pi_nl;

 %% save all images to tiff

 imwrite(imraw, 'unprocessed/imraw.tiff');
 imwrite(imconv, 'unprocessed/imconv.tiff')

 imwrite(im_vert_0pi_nl, 'unprocessed/im_vert_0pi_nl.tiff');
 imwrite(im_vert_0p4pi_nl, 'unprocessed/im_vert_0p4pi_nl.tiff');
 imwrite(im_vert_0p8pi_nl, 'unprocessed/im_vert_0p8pi_nl.tiff');
 imwrite(im_vert_1p2pi_nl, 'unprocessed/im_vert_1p2pi_nl.tiff');
 imwrite(im_vert_1p6pi_nl, 'unprocessed/im_vert_1p6pi_nl.tiff');

 imwrite(im_horz_0pi_nl, 'unprocessed/im_horz_0pi_nl.tiff');
 imwrite(im_horz_0p4pi_nl, 'unprocessed/im_horz_0p4pi_nl.tiff');
 imwrite(im_horz_0p8pi_nl, 'unprocessed/im_horz_0p8pi_nl.tiff');
 imwrite(im_horz_1p2pi_nl, 'unprocessed/im_horz_1p2pi_nl.tiff');
 imwrite(im_horz_1p6pi_nl, 'unprocessed/im_horz_1p6pi_nl.tiff');

 imwrite(im_diag_l45_0pi_nl, 'unprocessed/im_diag_l45_0pi_nl.tiff');
 imwrite(im_diag_l45_0p4pi_nl, 'unprocessed/im_diag_l45_0p4pi_nl.tiff');
 imwrite(im_diag_l45_0p8pi_nl, 'unprocessed/im_diag_l45_0p8pi_nl.tiff');
 imwrite(im_diag_l45_1p2pi_nl, 'unprocessed/im_diag_l45_1p2pi_nl.tiff');
 imwrite(im_diag_l45_1p6pi_nl, 'unprocessed/im_diag_l45_1p6pi_nl.tiff');

 imwrite(im_diag_r45_0pi_nl, 'unprocessed/im_diag_r45_0pi_nl.tiff');
 imwrite(im_diag_r45_0p4pi_nl, 'unprocessed/im_diag_r45_0p4pi_nl.tiff');
 imwrite(im_diag_r45_0p8pi_nl, 'unprocessed/im_diag_r45_0p8pi_nl.tiff');
 imwrite(im_diag_r45_1p2pi_nl, 'unprocessed/im_diag_r45_1p2pi_nl.tiff');
 imwrite(im_diag_r45_1p6pi_nl, 'unprocessed/im_diag_r45_1p6pi_nl.tiff');

 imwrite(psfraw, 'unprocessed/psfraw.tiff');

 imwrite(psfraw_vert_0pi_nl, 'unprocessed/psfraw_vert_0pi_nl.tiff');
 imwrite(psfraw_vert_0p4pi_nl, 'unprocessed/psfraw_vert_0p4pi_nl.tiff');
 imwrite(psfraw_vert_0p8pi_nl, 'unprocessed/psfraw_vert_0p8pi_nl.tiff');
 imwrite(psfraw_vert_1p2pi_nl, 'unprocessed/psfraw_vert_1p2pi_nl.tiff');
 imwrite(psfraw_vert_1p6pi_nl, 'unprocessed/psfraw_vert_1p6pi_nl.tiff');

 imwrite(psfraw_horz_0pi_nl, 'unprocessed/psfraw_horz_0pi_nl.tiff');
 imwrite(psfraw_horz_0p4pi_nl, 'unprocessed/psfraw_horz_0p4pi_nl.tiff');
 imwrite(psfraw_horz_0p8pi_nl, 'unprocessed/psfraw_horz_0p8pi_nl.tiff');
 imwrite(psfraw_horz_1p2pi_nl, 'unprocessed/psfraw_horz_1p2pi_nl.tiff');
 imwrite(psfraw_horz_1p6pi_nl, 'unprocessed/psfraw_horz_1p6pi_nl.tiff');

 imwrite(psfraw_diag_l45_0pi_nl, 'unprocessed/psfraw_diag_l45_0pi_nl.tiff');
 imwrite(psfraw_diag_l45_0p4pi_nl, 'unprocessed/psfraw_diag_l45_0p4pi_nl.tiff');
 imwrite(psfraw_diag_l45_0p8pi_nl, 'unprocessed/psfraw_diag_l45_0p8pi_nl.tiff');
 imwrite(psfraw_diag_l45_1p2pi_nl, 'unprocessed/psfraw_diag_l45_1p2pi_nl.tiff');
 imwrite(psfraw_diag_l45_1p6pi_nl, 'unprocessed/psfraw_diag_l45_1p6pi_nl.tiff');

 imwrite(psfraw_diag_r45_0pi_nl, 'unprocessed/psfraw_diag_r45_0pi_nl.tiff');
 imwrite(psfraw_diag_r45_0p4pi_nl, 'unprocessed/psfraw_diag_r45_0p4pi_nl.tiff');
 imwrite(psfraw_diag_r45_0p8pi_nl, 'unprocessed/psfraw_diag_r45_0p8pi_nl.tiff');
 imwrite(psfraw_diag_r45_1p2pi_nl, 'unprocessed/psfraw_diag_r45_1p2pi_nl.tiff');
 imwrite(psfraw_diag_r45_1p6pi_nl, 'unprocessed/psfraw_diag_r45_1p6pi_nl.tiff');

 imwrite(mat_vert_0pi, 'unprocessed/mat_vert_0pi.tiff');
 imwrite(mat_horz_0pi, 'unprocessed/mat_horz_0pi.tiff');
 imwrite(mat_diag_l45_0pi, 'unprocessed/mat_diag_l45_0pi.tiff');
 imwrite(mat_diag_r45_0pi, 'unprocessed/mat_diag_r45_0pi.tiff');
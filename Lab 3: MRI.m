% Name: Abby Yoon
% Pset 2, BE 128
% Date: 2/21/24

close all
Clc

% Name: Abby Yoon
% Pset 2, BE 128
% Date: 2/21/24
close all
clc

%% part 1. Reconstruct the four images (A1 through A4) from the raw k-space
A1 = load("A1.mat") % load the sinogram corresponding to an image
A2 = load("A2.mat") 
A3 = load("A3.mat") 
A4 = load("A4.mat") 
k_spaceA1 = A1.A1;
k_spaceA2 = A2.A2;
k_spaceA3 = A3.A3;
k_spaceA4 = A4.A4;

F1 = ifft2(k_spaceA1); % where k_spaceA1 is the raw k-space data
subplot(2,2,1), imagesc(F1), colormap(gray)
axis equal
axis off
title('full reconstruction of A1')

F2 = ifft2(k_spaceA2); 
subplot(2,2,2), imagesc(F2), colormap(gray)
axis equal
axis off
title('full reconstruction of A2')

F3 = ifft2(k_spaceA3); 
subplot(2,2,3), imagesc(F3), colormap(gray)
axis equal
axis off
title('full reconstruction of A3')

F4 = ifft2(k_spaceA4); 
subplot(2,2,4), imagesc(F4), colormap(gray)
axis equal
axis off
title('full reconstruction of A4')

%% number 2
% 2a - Plot your raw k-space data
b1 = abs(fftshift(k_spaceA1));
figure, imagesc(b1), colormap(gray)
axis equal
axis off
title('raw k-space')

% 2b - Compressed vs uncompressed k-space data
b2 = 20*log10(b1/(max(max(b1))));
figure, imagesc(b2), colormap(gray)
axis equal
axis off
title('raw k-space compressed')

% 2c - Reconstruct the image using only the center quarter of the k-space
% (128 by 128 points)
B1 = fftshift(k_spaceA1);
filt = zeros(size(B1)); % makes all signals zeroes
[ph,fr] = size(B1); % phase encoding in x-direction and frequency encoding in y-direction
x = ph/2+0; % this is the center of the phase direction of k-space
y = fr/2+0; % this is the center of the frequency direction of k-space
w = 128;
h = 128;
filt(round(y-h/2):round(y+h/2),round(x-w/2):round(x+w/2)) = 1; % takes all the zeros in this range and converts them to ones (throwing away the outer points)
conv = filt.*B1; % multiply the filter by the original k-space (convolution)
C = fftshift(conv); % uncenter it to prepare for 2-D inverse Fourier Transform
D = ifft2(C, 'symmetric');
figure, imagesc(abs(D)), colormap(gray)
axis off
axis equal

% 2d -  Reconstruct the image using only the center quarter of the k-space
% (64 by 64 points)
[ph,fr] = size(B1); % phase encoding in x-direction and frequency encoding in y-direction
x = ph/2+0; % this is the center of the phase direction of k-space
y = fr/2+0; % this is the center of the frequency direction of k-space
w = 64;
h = 64;
filt(round(y-h/2):round(y+h/2),round(x-w/2):round(x+w/2)) = 1; % takes all the zeros in this range and converts them to ones (throwing away the outer points)
conv = filt.*B1; % multiply the filter by the original k-space (convolution)
C = fftshift(conv); % uncenter it to prepare for 2-D inverse Fourier Transform
D = ifft2(C, 'symmetric');
figure, imagesc(abs(D)), colormap(gray)
axis off
axis equal

% 2e -  Reconstruct the image using only the outer part of the k-space
% (the rest of the matrix after you have thrown away 128 by 128 points)
B1 = fftshift(k_spaceA1);
filt = ones(size(B1)); % makes all the signals ones
[ph,fr] = size(B1); % phase and frequency encoding
x = ph/2+0; % this is the center of the phase direction of k-space
y = fr/2+0; % this is the center of the frequency direction of k-space
w = 128;
h = 128;
filt(round(y-h/2):round(y+h/2),round(x-w/2):round(x+w/2)) = 0; % takes all the ones in this range and converts them to zeroes (throwing away the center points)
conv = filt.*B1; % multiply the filter by the original k-space (convolution)
C = fftshift(conv); % uncenter it to prepare for 2-D inverse Fourier Transform
D = ifft2(C, 'symmetric');
figure, imagesc(abs(D)), colormap(gray)
axis off
axis equal

% 2f - Remove half of the k-space points and reconstruct the image from your
% downsampled k-space.
B1 = fftshift(k_spaceA1); % centered
DD = downsample(B1,2); % takes every other point in the x-direction
D3 = DD'; % transpose the matrix to get the y-values as x-values
D4 = downsample(D3,2); % takes every other point in the y-direction
D4 = D4'; % transpose the matrix again
C = fftshift(D4); % uncenter it to prepare for inverse FT
D = ifft2(C, 'symmetric');
figure, imagesc(abs(D)), colormap(gray)
axis off
axis equal

%% number 3
% 3A Reconstruct the image using only a single, off-center, point of the k-space
% % create a simple filter of ones and zeros with a single k space element
A2 = fftshift(k_spaceA1); %centered
B = zeros(size(A2));
[ph,fr] = size(A2);
x = ph/2+10;   % this is the center of the phase direction of k-space
%offset by one (change the amount added to shift further from the center
%of k-space)
y = fr/2+5; % this is the center of the frequency direction of k-space
%offset by 5 (change the amount added to shift further from the center
%of k-space)
B(x,y)= 1;
filt = B;
%figure
%subplot(221)
%imagesc(filt), axis off, axis equal
conv = filt.*A2;
C = fftshift(conv); %uncenter it to prepare for inverse FT
D = ifft2(C, 'symmetric');
figure
% subplot(221) % plot the image from the filtered k-space
imagesc(abs(D)), axis off, axis equal
colormap(gray)


% 3B Reconstruct the image using only a 16x 64 portion of the k-space by creating a filter of ones
% and zeros (ones inside the 16x64 space and zeros everywhere else)
B1 = fftshift(k_spaceA1);
[ph,fr] = size(B1); % phase encoding in x-direction and frequency encoding in y-direction
x = ph/2+0; % this is the center of the phase direction of k-space
y = fr/2+0; % this is the center of the frequency direction of k-space
w = 16;
h = 64;
filt(round(y-h/2):round(y+h/2),round(x-w/2):round(x+w/2)) = 1; % takes all the zeros in this range and converts them to ones (throwing away the outer points)
conv = filt.*B1; % multiply the filter by the original k-space (convolution)
C = fftshift(conv); % uncenter it to prepare for 2-D inverse Fourier Transform
D = ifft2(C, 'symmetric');
figure, imagesc(abs(D)), colormap(gray)
axis off
axis equal


% 3C Reconstruct the image using only the pixels of the k-space that are above a threshold. I
% suggest using the top 10% or 25% brightest pixels
B1 = fftshift(k_spaceA1);
filt = zeros(size(B1));
k_bs = abs(B1);
threshold = .75*(max(max(k_bs)));
above_threshold = k_bs > threshold;
conv = above_threshold.*B1
C = fftshift(conv); % uncenter it to prepare for 2-D inverse Fourier Transform
D = ifft2(C, 'symmetric');
figure
imagesc(abs(D)), colormap(gray), axis equal, axis off


% 3D Create a filter that only keeps the high frequencies (just like in part 1e) of k-space.
% Adjust your filter until you have an image that mostly looks like
% edges.
B1 = fftshift(k_spaceA1);
filt = ones(size(B1)); % makes all the signals ones
[ph,fr] = size(B1); % phase and frequency encoding
x = ph/2+0; % this is the center of the phase direction of k-space
y = fr/2+0; % this is the center of the frequency direction of k-space
w = 205;
h = 205;
filt(round(y-h/2):round(y+h/2),round(x-w/2):round(x+w/2)) = 0; % takes all the ones in this range and converts them to zeroes (throwing away the center points)
conv = filt.*B1; % multiply the filter by the original k-space (convolution)
C = fftshift(conv); % uncenter it to prepare for 2-D inverse Fourier Transform
D = ifft2(C, 'symmetric');
figure, imagesc(abs(D)), colormap(gray)
axis off
axis equal


%% BONUS exercise #2: Edge detection from sharper frequencies
% Sobel
BW1 = edge(F1); % F1 was defined many parts ago - F1 = ifft2(k_spaceA1)
subplot(2,3,1)
imshow(BW1)
title('Sobel')
% Prewitt
BW2 = edge(F1,'Prewitt');
subplot(2,3,2)
imshow(BW2)
title('Prewitt')
% Roberts
BW3 = edge(F1,'Roberts');
subplot(2,3,3)
imshow(BW3);
title('Roberts')
% log
BW4 = edge(F1,'log');
subplot(2,3,4)
imshow(BW4);
title('Log')
% zerocross
BW5 = edge(F1,'zerocross');
subplot(2,3,5)
imshow(BW5);
title('Zero-Cross')
% canny
BW6 = edge(F1,'Canny');
subplot(2,3,6)
imshow(BW6);
title('Canny')


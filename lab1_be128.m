% M-file for focusing the raw RF ultrasound data
% Functions called      : calc_apod_sonos.m
%                         log_comp.m
% Date of Original Creation      : 5/8/06
% Original Authors: Drake Guenther, Linsey Phillips, Priya Raghunathan
% Date of last edit:    1/31/13

%% Reconstruction for cysts

clear all
load rf_cysts          %loads rf data file to be used
d = rf_cysts;      % creates a variable with the values of the data file    
clear rf_*

% Known input parameters
fs = 39.27e6;           % sampling frequency
T = 1/fs;               % sampling period
ele_pitch =  135e-6;    % distance (pitch) between elements
c = 1540;               % m/s this is assumed for the phantom
num_ele = 128;          % number of elements in the transducer
lat_dx = ele_pitch*[-64:63]+ele_pitch/2;    %lateral dimensions
[d_x d_y] = size(d);
ax_t = [0:d_x-1]*T;                        
t_offset = 174*T;                           %deadtime # of samples
ax_t = ax_t+t_offset;    %axial dimensions with deadtime taken into account
ax_d = ax_t*c/2;   %divide by two for pulse-echo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Interpolation of data and new vectors %%%%%%%%%%%%
interp_fac = 4;         %interpolates the raw data
no_lines = 128*interp_fac-1;    %creates new lines of data
beam_spacing = ele_pitch/interp_fac;
line_dx = beam_spacing*[-no_lines:no_lines];
[X,Z] = meshgrid(line_dx,ax_d);
del_tx = Z/c;                    %delay times on transmit
del_rx = realsqrt(X.^2+Z.^2)/c;  %delay times on receive
foc_time = del_tx+del_rx;        %total delay time for each element
               
%generates an FIR filter to remove noise from rf_data spectrum
ord = 300;
f = fir1(ord,[(2)/(39/2) (9)/(39/2)]); 
%actual filtering done by convolution
for n = 1:128
   df(:,n) = conv(d(:,n),f);
end
df = df(151:end-150,:);  
clear X Z d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CALCULATE APODIZATION PROFILE %%%%%%%%%%%%
fn = 0.75;  % f-number used on rx
w=calc_apod_sonos(fn,ax_d,interp_fac);
%w = ones(size(foc_time,1),size(foc_time),2);   % Rect window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redefine the raw echo data in terms of splines
for j=1:num_ele
   eval(['dpp' int2str(j) '=spline(ax_t,df(:,j));']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pack
% define the matrix that will hold the delayed echo data
rf=zeros(length(ax_d),num_ele,no_lines,'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  FOCUS DATA BY EVALUATING SPLINES AT FOC TIMES %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
as=no_lines+1;
for j=1:num_ele
   for k=1:no_lines;
       eval(['rf(:,j,k)=ppval(dpp' int2str(j) ...
           ', foc_time(:,as+k-1)).*w(:,as+k-1);']);
   end
   as=as-interp_fac;
   j/num_ele*100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beamform (sum) for alines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alines = zeros(2100,no_lines);
lat_dist = line_dx/2;
for n = 1:no_lines
   alines(:,n) = sum(rf(1:2100,:,n),2);
end
alines = alines(200:2000,:);
figure, imagesc(lat_dist*1000,ax_d*1000,alines), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Time Delayed Image')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Envelope detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
env = abs(hilbert(alines(:,:)));
figure
imagesc(lat_dist*1e3,ax_d*1e3,env), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Enveloped Image')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log Compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logenv = log_comp(env,40);
figure
imagesc(lat_dist*1e3,ax_d*1e3, logenv), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Log Compressed')
%% Plot original image
figure
imagesc(lat_dist*1e3,ax_d*1e3,d), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Original Image')

%% Reconstruction for legions

     % for rf_lesions
clear all
load rf_lesions          %loads rf data file to be used
d = rf_lesions;      % creates a variable with the values of the data file    
clear rf_*

% Known input parameters
fs = 39.27e6;           % sampling frequency
T = 1/fs;               % sampling period
ele_pitch =  135e-6;    % distance (pitch) between elements
c = 1540;               % m/s this is assumed for the phantom
num_ele = 128;          % number of elements in the transducer
lat_dx = ele_pitch*[-64:63]+ele_pitch/2;    %lateral dimensions
[d_x d_y] = size(d);
ax_t = [0:d_x-1]*T;                        
t_offset = 174*T;                           %deadtime # of samples
ax_t = ax_t+t_offset;    %axial dimensions with deadtime taken into account
ax_d = ax_t*c/2;   %divide by two for pulse-echo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Interpolation of data and new vectors %%%%%%%%%%%%
interp_fac = 4;         %interpolates the raw data
no_lines = 128*interp_fac-1;    %creates new lines of data
beam_spacing = ele_pitch/interp_fac;
line_dx = beam_spacing*[-no_lines:no_lines];
[X,Z] = meshgrid(line_dx,ax_d);
del_tx = Z/c;                    %delay times on transmit
del_rx = realsqrt(X.^2+Z.^2)/c;  %delay times on receive
foc_time = del_tx+del_rx;        %total delay time for each element
               
%generates an FIR filter to remove noise from rf_data spectrum
ord = 300;
f = fir1(ord,[(2)/(39/2) (9)/(39/2)]); 
%actual filtering done by convolution
for n = 1:128
   df(:,n) = conv(d(:,n),f);
end
df = df(151:end-150,:);  
clear X Z d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CALCULATE APODIZATION PROFILE %%%%%%%%%%%%
fn = 0.75;  % f-number used on rx
w=calc_apod_sonos(fn,ax_d,interp_fac);
%w = ones(size(foc_time,1),size(foc_time),2);   % Rect window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redefine the raw echo data in terms of splines
for j=1:num_ele
   eval(['dpp' int2str(j) '=spline(ax_t,df(:,j));']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pack
% define the matrix that will hold the delayed echo data
rf=zeros(length(ax_d),num_ele,no_lines,'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  FOCUS DATA BY EVALUATING SPLINES AT FOC TIMES %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
as=no_lines+1;
for j=1:num_ele
   for k=1:no_lines;
       eval(['rf(:,j,k)=ppval(dpp' int2str(j) ...
           ', foc_time(:,as+k-1)).*w(:,as+k-1);']);
   end
   as=as-interp_fac;
   j/num_ele*100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beamform (sum) for alines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alines = zeros(2100,no_lines);
lat_dist = line_dx/2;
for n = 1:no_lines
   alines(:,n) = sum(rf(1:2100,:,n),2);
end
alines = alines(200:2000,:);
figure, imagesc(lat_dist*1000,ax_d*1000,alines), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Time Delayed Image')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Envelope detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
env = abs(hilbert(alines(:,:)));
figure
imagesc(lat_dist*1e3,ax_d*1e3,env), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Enveloped Image')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log Compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logenv = log_comp(env,40);
figure
imagesc(lat_dist*1e3,ax_d*1e3, logenv), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Log Compressed')
%% Reconstruction for wires
clear all
load rf_wires          %loads rf data file to be used
d = rf_wires;      % creates a variable with the values of the data file    
clear rf_*
% Known input parameters
fs = 39.27e6;           % sampling frequency
T = 1/fs;               % sampling period
ele_pitch =  135e-6;    % distance (pitch) between elements
c = 1540;               % m/s this is assumed for the phantom
num_ele = 128;          % number of elements in the transducer
lat_dx = ele_pitch*[-64:63]+ele_pitch/2;    %lateral dimensions
[d_x d_y] = size(d);
ax_t = [0:d_x-1]*T;                        
t_offset = 174*T;                           %deadtime # of samples
ax_t = ax_t+t_offset;    %axial dimensions with deadtime taken into account
ax_d = ax_t*c/2;   %divide by two for pulse-echo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Interpolation of data and new vectors %%%%%%%%%%%%
interp_fac = 4;         %interpolates the raw data
no_lines = 128*interp_fac-1;    %creates new lines of data
beam_spacing = ele_pitch/interp_fac;
line_dx = beam_spacing*[-no_lines:no_lines];
[X,Z] = meshgrid(line_dx,ax_d);
del_tx = Z/c;                    %delay times on transmit
del_rx = realsqrt(X.^2+Z.^2)/c;  %delay times on receive
foc_time = del_tx+del_rx;        %total delay time for each element
               
%generates an FIR filter to remove noise from rf_data spectrum
ord = 300;
f = fir1(ord,[(2)/(39/2) (9)/(39/2)]); 
%actual filtering done by convolution
for n = 1:128
   df(:,n) = conv(d(:,n),f);
end
df = df(151:end-150,:);  
clear X Z d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CALCULATE APODIZATION PROFILE %%%%%%%%%%%%
fn = 0.75;  % f-number used on rx
w=calc_apod_sonos(fn,ax_d,interp_fac);
%w = ones(size(foc_time,1),size(foc_time),2);   % Rect window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Redefine the raw echo data in terms of splines
for j=1:num_ele
   eval(['dpp' int2str(j) '=spline(ax_t,df(:,j));']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pack
% define the matrix that will hold the delayed echo data
rf=zeros(length(ax_d),num_ele,no_lines,'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  FOCUS DATA BY EVALUATING SPLINES AT FOC TIMES %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
as=no_lines+1;
for j=1:num_ele
   for k=1:no_lines;
       eval(['rf(:,j,k)=ppval(dpp' int2str(j) ...
           ', foc_time(:,as+k-1)).*w(:,as+k-1);']);
   end
   as=as-interp_fac;
   j/num_ele*100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beamform (sum) for alines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alines = zeros(2100,no_lines);
lat_dist = line_dx/2;
for n = 1:no_lines
   alines(:,n) = sum(rf(1:2100,:,n),2);
end
alines = alines(200:2000,:);
figure, imagesc(lat_dist*1000,ax_d*1000,alines), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Time Delayed Image')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Envelope detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
env = abs(hilbert(alines(:,:)));
figure
imagesc(lat_dist*1e3,ax_d*1e3,env), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Enveloped Image')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log Compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logenv = log_comp(env,40);
figure
imagesc(lat_dist*1e3,ax_d*1e3, logenv), colormap(gray), axis equal tight
xlabel('Lateral distance (mm)')
ylabel('Axial distance (mm)')
title('Log Compressed')

%% PART 1 POST-LAB CALCULATIONS

% Calculate fundamental frequency of ultrasound signals in 'rf_wires'
   
   % load rf_wires data
   d = load('rf_wires.mat');  
   fs = 39.27e6;           % sampling frequency
   N0 = 39.27e6;
   data  = d.rf_wires;
   data = reshape(data, 1, []);
   X = fftshift(abs(fft(data,N0))); % shifting is needed to have our data reflect the symmetrical nature of our data between its positive and negative values of the frequency domain
   % make x-axis for plotting the signal
   F0 = ([-N0/2:N0/2-1]/N0)*fs;
   % plot fourier
   figure;
   plot(F0, X)
   title("FT of Signal X (Shifted)")
   ylabel("Amplitude [A.u.]")
   xlabel("Frequency [Hz]")

 % Calculate the theoretical lateral resolution of linear array transducer
   % lateral resolution at 6.5 cm depth: F-Number * Wavelength =  λ*L/D
    width = 4; % in 'cm'(D)
    focal_length = 6.5; % (L)
    c = 1540; % speed of sound (c), don't need to convert m/s to cm/s
    f = 800000; % frequency (f) in Hz, this is = 8 MHz 
    wavelength = c/f;
    lateral_resolution = wavelength * focal_length / width
    % lateral resolution at 3 cm depth: F-Number * Wavelength =  λ*L/D
    width = 4; % in 'cm'(D)
    focal_length = 3; % (L)
    c = 1540; % speed of sound (c), don't need to convert m/s to cm/s
    f = 800000; % frequency (f) in Hz, this is = 8 MHz 
    wavelength = c/f;
    lateral_resolution = wavelength * focal_length / width

% Calculate peak negative pressure in Q12
    mechanical_index = 0.9; % in MPa/ MHz^-1/2
    frequency = 10; % in MHz
    pnp = mechanical_index * (frequency^(-1/2))


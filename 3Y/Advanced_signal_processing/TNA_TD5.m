% TD5 TNS Avancé DU
% TD Machine: Design of filter-banks
%%
close all ;
clear;
clc ;
format short g ;

%% Design a tree-structure 8-channel filter bank.
% Design the PR Orthogonal 2-ch filter bank

% Corrected Design parameters
N = 33;                   % Filter length (must be odd)
wp = 0.43;                % Normalized lowpass passband edge frequency (0 < wp < 0.5)

% Design the two-channel orthogonal filter bank using firpr2chfb
[h0, h1, g0, g1] = firpr2chfb(N, wp);

% Plot the impulse responses of analysis and synthesis filters
figure;
subplot(2,2,1);
stem(h0);
title('Impulse Response h_0[n] (Analysis Lowpass)');
xlabel('n');
ylabel('h_0[n]');

subplot(2,2,2);
stem(h1);
title('Impulse Response h_1[n] (Analysis Highpass)');
xlabel('n');
ylabel('h_1[n]');

subplot(2,2,3);
stem(g0);
title('Impulse Response g_0[n] (Synthesis Lowpass)');
xlabel('n');
ylabel('g_0[n]');

subplot(2,2,4);
stem(g1);
title('Impulse Response g_1[n] (Synthesis Highpass)');
xlabel('n');
ylabel('g_1[n]');

% Compute and plot poles and zeros of transfer functions H0(z), H1(z), G0(z), G1(z)
figure;
subplot(2,2,1);
zplane(h0, 1);
title('Poles and Zeros of H_0(z)');

subplot(2,2,2);
zplane(h1, 1);
title('Poles and Zeros of H_1(z)');

subplot(2,2,3);
zplane(g0, 1);
title('Poles and Zeros of G_0(z)');

subplot(2,2,4);
zplane(g1, 1);
title('Poles and Zeros of G_1(z)');

% Compute and plot magnitude and group delay responses of analysis filters H0(z) and H1(z)
[H0, w] = freqz(h0, 1, 1024);
[H1, w] = freqz(h1, 1, 1024);
grp_delay_h0 = grpdelay(h0, 1, 1024);
grp_delay_h1 = grpdelay(h1, 1, 1024);

figure;
subplot(2,2,1);
plot(w/pi, abs(H0));
title('Magnitude Response of H_0(z)');
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('|H_0(e^{j\omega})|');

subplot(2,2,2);
plot(w/pi, abs(H1));
title('Magnitude Response of H_1(z)');
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('|H_1(e^{j\omega})|');

subplot(2,2,3);
plot(w/pi, grp_delay_h0);
title('Group Delay of H_0(z)');
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Group Delay (samples)');

subplot(2,2,4);
plot(w/pi, grp_delay_h1);
title('Group Delay of H_1(z)');
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Group Delay (samples)');

% Compute and plot the impulse response of the distortion transfer function t[n]
t = conv(h0, g0) + conv(h1, g1);

figure;
stem(t);
title('Impulse Response of Distortion Transfer Function t[n]');
xlabel('n');
ylabel('t[n]');

%%

% Filter length and normalized passband edge frequency
N = 33;  % Filter length (must be odd)
wp = 0.43;  % Normalized lowpass passband edge frequency (0 < wp < 0.5)

% Design the two-channel orthogonal filter bank using firpr2chfb
[h0, h1, g0, g1] = firpr2chfb(N, wp);

% Define the number of channels
num_channels = 8;

% Prepare to store the eight analysis filters
filters = cell(num_channels, 1);

% Perform the tree-structured decomposition for the analysis filter bank
% Level 1: Apply two-channel filter bank on original signal
filters{1} = h0;  % Lowpass branch
filters{2} = h1;  % Highpass branch

% Level 2: Apply two-channel filter bank on each of the above branches
filters{3} = conv(h0, h0);  % Lowpass filter applied to Lowpass branch
filters{4} = conv(h0, h1);  % Highpass filter applied to Lowpass branch
filters{5} = conv(h1, h0);  % Lowpass filter applied to Highpass branch
filters{6} = conv(h1, h1);  % Highpass filter applied to Highpass branch

% Level 3: Apply two-channel filter bank on each of the level 2 branches
filters{7} = conv(filters{3}, h1);  % Further decompose branch 3 (low-low)
filters{8} = conv(filters{4}, h1);  % Further decompose branch 4 (low-high)

% Plot magnitude responses of the resulting eight analysis filters
figure;
for k = 1:num_channels
    [Hk, w] = freqz(filters{k}, 1, 1024);  % Compute frequency response
    subplot(4, 2, k);  % Plot 4x2 grid
    plot(w/pi, abs(Hk));  % Plot magnitude response
    title(['Magnitude Response of Analysis Filter ', num2str(k)]);
    xlabel('Normalized Frequency (\omega/\pi)');
    ylabel(['|H', num2str(k), '(e^{j\omega})|']);
end

% Verify the magnitude-preserving property of the overall bank
total_response = zeros(size(w));
for k = 1:num_channels
    [Hk, w] = freqz(filters{k}, 1, 1024);
    total_response = total_response + abs(Hk).^2;  % Sum of magnitudes squared
end

% Plot the total response
figure;
plot(w/pi, total_response);
title('Total Magnitude Response of the Filter Bank');
xlabel('Normalized Frequency (\omega/\pi)');
ylabel('Sum of Squared Magnitudes');

% Check if the total response is approximately 1 across all frequencies
if all(abs(total_response - 1) < 1e-5)
    disp('The filter bank satisfies the magnitude-preserving property.');
else
    disp('The filter bank does not preserve magnitude.');
end

%%
% Create the rectangular signal: 200 samples of 1 followed by 1000 samples of 0
rect_signal = [ones(1, 200), zeros(1, 1000)];

% Filter length and normalized passband edge frequency
N = 33;  % Filter length (must be odd)
wp = 0.43;  % Normalized lowpass passband edge frequency (0 < wp < 0.5)

% Design the two-channel orthogonal filter bank using firpr2chfb
[h0, h1, g0, g1] = firpr2chfb(N, wp);

% Prepare to store the eight analysis filters
filters = cell(8, 1);

% Perform the tree-structured decomposition for the analysis filter bank
filters{1} = h0;  % Lowpass branch
filters{2} = h1;  % Highpass branch
filters{3} = conv(h0, h0);  % Lowpass filter applied to Lowpass branch
filters{4} = conv(h0, h1);  % Highpass filter applied to Lowpass branch
filters{5} = conv(h1, h0);  % Lowpass filter applied to Highpass branch
filters{6} = conv(h1, h1);  % Highpass filter applied to Highpass branch
filters{7} = conv(filters{3}, h1);  % Further decompose branch 3 (low-low)
filters{8} = conv(filters{4}, h1);  % Further decompose branch 4 (low-high)

% Decompose the rectangle signal using the eight filters
channel_outputs = cell(8, 1);
for k = 1:8
    channel_outputs{k} = filter(filters{k}, 1, rect_signal);
end

% Plot the time responses of all eight channels
figure;
for k = 1:8
    subplot(4, 2, k);  % 4x2 grid for plotting
    plot(channel_outputs{k});
    title(['Time Response of Channel ', num2str(k)]);
    xlabel('Sample index');
    ylabel('Amplitude');
end

% Check the SBC (Subband Coding) Impulse Response

% Create an impulse signal
impulse_signal = [1, zeros(1, 1199)];  % Impulse of length 1200 samples

% Pass the impulse signal through the filter bank
impulse_outputs = cell(8, 1);
for k = 1:8
    impulse_outputs{k} = filter(filters{k}, 1, impulse_signal);
end

% Plot the impulse responses of all eight channels
figure;
for k = 1:8
    subplot(4, 2, k);  % 4x2 grid for plotting
    stem(impulse_outputs{k});
    title(['Impulse Response of Channel ', num2str(k)]);
    xlabel('Sample index');
    ylabel('Amplitude');
end



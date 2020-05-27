%%% Yuval Epstain Ofek 
%%% Sound and Space - Problem Set 6

%%
clear all; close all; clc;
%Some parameters to generate noise
N = 1000000;
fs = 10000;
Pref = 2e-5;

%%% a.
%Generating noise
p = pinknoise(N,fs);
w = whitenoise(N,fs);


%Finding RMS and ensuring both signals have the same RMS (chose the RMS of
%the white noise arbitrarily
RMSp = sqrt(mean(p.^2))
RMSw = sqrt(mean(w.^2))
w = w/RMSw*RMSp;

%checking that they are in fact the same
RMSp = sqrt(mean(p.^2))
RMSw = sqrt(mean(w.^2))

%plotting the generated time series
figure 
subplot(2,1,1)
plot((0:N-1)/fs, real(w), 'LineWidth', 1.5)
title('White Noise','FontSize', 16)
xlabel('Time [s]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Noise', 'FontSize', 16)
xlim([0,(N-1)/fs])
subplot(2,1,2)
plot((0:N-1)/fs, real(p), 'LineWidth', 1.5)
title('Pink Noise','FontSize', 16)
xlabel('Time [s]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Noise', 'FontSize', 16)
xlim([0,(N-1)/fs])
sgtitle('Time series of the noise generated','FontWeight', 'bold', 'FontSize', 20)

%% b. 
close all; clc;
%Since both signals are normalized to have the same RMS value, they should
%also have the same overall SPL (as SPL is simply the 20log10(RMS/P_ref))
SPLdBw = 20*log10(RMSw/Pref)
SPLdBp = 20*log10(RMSp/Pref)

%% c.
clc

%Getting our noise to the frequency domain
[W,f] = my_fft(w,N,fs);
[P,~] = my_fft(p,N,fs); %since fs and N are the same for p & w, the frequencies will be too.

%Frequencies to get filters
f1 = 12194.22;
f2 = 20.598997;
f3 = 107.65265;
f4 = 737.86223;
%anon function to convert from Hz to rad/s
f2w = @(f) 2*pi*f;
%converting the frequencies
w1 = f2w(f1);
w2 = f2w(f2);
w3 = f2w(f3);
w4 = f2w(f4);
wvect = f2w(f); 

%Generating filter H_a & H_c at the frequencies we want
Hc = 1.0072*(1j*wvect).^2*w1^2./((w1+1j*wvect).^2.*(w2+1j*wvect).^2);
Ha = Hc.*(1.25*(1j*wvect).^2./((w3+1j*wvect).*(w4+1j*wvect)));


%Computing the weighted noise 
Wc = W.*Hc;
Wa = W.*Ha;
Pc = P.*Hc;
Pa = P.*Ha;
w_c = my_ifft(Wc,N,fs);
p_c = my_ifft(Pc,N,fs);

w_a = my_ifft(Wa,N,fs);
p_a = my_ifft(Pa,N,fs);

%plotting
%time series
figure 
subplot(2,2,1)
plot((0:N-1)/fs, real(w_c), 'LineWidth', 1.5)
title('White Noise C-Weighting','FontSize', 16)
xlabel('Time [s]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Noise', 'FontSize', 16)
xlim([0,(N-1)/fs])

subplot(2,2,3)
plot((0:N-1)/fs, real(p_c), 'LineWidth', 1.5)
title('Pink Noise C-Weighting','FontSize', 16)
xlabel('Time [s]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Noise', 'FontSize', 16)
xlim([0,(N-1)/fs])

subplot(2,2,2)
plot((0:N-1)/fs, real(w_a), 'LineWidth', 1.5)
title('White Noise A-Weighting','FontSize', 16)
xlabel('Time [s]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Noise', 'FontSize', 16)
xlim([0,(N-1)/fs])

subplot(2,2,4)
plot((0:N-1)/fs, real(p_a), 'LineWidth', 1.5)
title('Pink Noise A-Weighting','FontSize', 16)
xlabel('Time [s]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Noise', 'FontSize', 16)
xlim([0,(N-1)/fs])
sgtitle('Time series of the weighted noise','FontWeight', 'bold', 'FontSize', 20)


%%% SPL calculation
SPLWcdB = 10*log(mean(w_c.^2)/Pref^2)
SPLWadB = 10*log(mean(w_a.^2)/Pref^2)
SPLPcdB = 10*log(mean(p_c.^2)/Pref^2)
SPLPadB = 10*log(mean(p_a.^2)/Pref^2)

%% d. 
clc; close all;
%creating a vector of the cuttoff frequencies
exp = 0:5;
f_cut = 125*2.^exp;
% to normalized frequency
w_cut = f_cut*2/fs;

%initialize 
Ba = zeros(5,7);
Bb = zeros(5,7);
w_oct = zeros(5, N);
W_OCT = zeros(5, N);
p_oct = zeros(5, N);
P_OCT = zeros(5, N);

%Get filter coefficients
Filters = 1:size(f_cut, 2)-1;
for i = Filters
    [Bb(i,:),Ba(i,:)] = butter(3, [w_cut(i),w_cut(i+1)], 'bandpass');
    w_oct(i,:) = filter(Bb(i,:),Ba(i,:),w);
    [W_OCT(i,:),f] = my_fft(w_oct(i,:),N,fs);  %freq resp
    p_oct(i,:) = filter(Bb(i,:),Ba(i,:),p);
    P_OCT(i,:) = my_fft(p_oct(i,:),N,fs);  %freq resp
end

%%% e.
%initialize 
wa_oct = zeros(5, N);
pa_oct = zeros(5, N);
wc_oct = zeros(5, N);
pc_oct = zeros(5, N);

for i = Filters
    wa_oct(i,:) = filter(Bb(i,:),Ba(i,:),w_a);
    wc_oct(i,:) = filter(Bb(i,:),Ba(i,:),w_c);
    pa_oct(i,:) = filter(Bb(i,:),Ba(i,:),p_a);
    pc_oct(i,:) = filter(Bb(i,:),Ba(i,:),p_c);
end
% An alternate way to do this, I like the first one more because I felt
% like it was more informative. 
%{
for i = Filters
    wa_oct(i,:) = my_ifft(W_OCT(i,:).*Ha,N,fs);
    wc_oct(i,:) = my_ifft(W_OCT(i,:).*Hc,N,fs);
    pa_oct(i,:) = my_ifft(P_OCT(i,:).*Ha,N,fs);
    pc_oct(i,:) = my_ifft(P_OCT(i,:).*Hc,N,fs);
end
%}
SPL_w = 10*log10(mean(w_oct.^2,2)/Pref^2);
SPL_p = 10*log10(mean(p_oct.^2,2)/Pref^2);

SPL_wa = 10*log10(mean(wa_oct.^2,2)/Pref^2);
SPL_wc = 10*log10(mean(wc_oct.^2,2)/Pref^2);
SPL_pa = 10*log10(mean(pa_oct.^2,2)/Pref^2);
SPL_pc = 10*log10(mean(pc_oct.^2,2)/Pref^2);

figure 
subplot(2,1,1)
semilogx(Filters, real(real(SPL_w)), '-o',Filters, real(real(SPL_wa)), '-o', ...
    Filters, real(real(SPL_wc)), '-o','LineWidth', 1.5)
title('White Noise','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('SPL [dB re: 20\muPa]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
xticks(Filters)
xticklabels({'125 - 250', '250 - 500', '500 - 1,000', '1,000 - 2,000', '2,000 - 4,000'})
legend('Overall SPL', 'A-Weighted', 'C-Weighted', 'FontSize', 16)

subplot(2,1,2)
semilogx(Filters, real(real(SPL_p)), '-o',Filters, real(real(SPL_pa)), '-o', ...
    Filters, real(real(SPL_pc)), '-o','LineWidth', 1.5)
title('Pink Noise','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('SPL [dB re: 20\muPa]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
xticks(Filters)
xticklabels({'125 - 250', '250 - 500', '500 - 1,000', '1,000 - 2,000', '2,000 - 4,000'})
legend('Overall SPL', 'A-Weighted', 'C-Weighted', 'FontSize', 16)
sgtitle('Level vs. Octave plot','FontWeight', 'bold', 'FontSize', 20)

% ii)
%Freq Response
WcdB = 20*log(Wc/Pref);
WadB = 20*log(Wa/Pref);
PcdB = 20*log(Pc/Pref);
PadB = 20*log(Pa/Pref);

figure 
subplot(2,1,1)
sgtitle('Frequency response of the noise','FontWeight', 'bold', 'FontSize', 20)
semilogx(-f(1:size(f,2)/2),real(20*log10(W(1:size(f,2)/2)/Pref)), ...
    -f(1:size(f,2)/2), real(WcdB(1:size(f,2)/2)),-f(1:size(f,2)/2), ... 
    real(WadB(1:size(f,2)/2)), 'LineWidth', 1.5)
title('White Noise','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [dB re: 20\muPa]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
xlim([0,fs/2])
legend('Non-Weighted', 'C-Weighting','A-Weighting', 'FontSize', 16)


subplot(2,1,2)
semilogx(-f(1:size(f,2)/2),real(20*log10(P(1:size(f,2)/2)/Pref)), ...
    -f(1:size(f,2)/2), real(PcdB(1:size(f,2)/2)),-f(1:size(f,2)/2), ...
    real(PadB(1:size(f,2)/2)), 'LineWidth', 1.5)
title('Pink Noise','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [dB re: 20\muPa]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Non-Weighted', 'C-Weighting','A-Weighting', 'FontSize', 16)
xlim([0,fs/2])

%%% Plotting again, but focused on the frequencies we want:
figure 
subplot(2,1,1)
sgtitle('Frequency response of the noise','FontWeight', 'bold', 'FontSize', 20)
semilogx(-f(1:size(f,2)/2),real(20*log10(W(1:size(f,2)/2)/Pref)), ...
    -f(1:size(f,2)/2), real(WcdB(1:size(f,2)/2)),-f(1:size(f,2)/2), ... 
    real(WadB(1:size(f,2)/2)), 'LineWidth', 1.5)
title('White Noise','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [dB re: 20\muPa]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
xlim([125,4000])
legend('Non-Weighted', 'C-Weighting','A-Weighting', 'FontSize', 16)
for i = 1:size(f_cut,2)
    x = xline(f_cut(i), '--', {[num2str(f_cut(i)),' Hz']}, 'LineWidth', 2);
    set(get(get(x,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

subplot(2,1,2)
semilogx(-f(1:size(f,2)/2),real(20*log10(P(1:size(f,2)/2)/Pref)), ...
    -f(1:size(f,2)/2), real(PcdB(1:size(f,2)/2)),-f(1:size(f,2)/2), ...
    real(PadB(1:size(f,2)/2)), 'LineWidth', 1.5)
title('Pink Noise','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [dB re: 20\muPa]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Non-Weighted', 'C-Weighting','A-Weighting', 'FontSize', 16)
xlim([125,4000])
for i = 1:size(f_cut,2)
    x = xline(f_cut(i), '--', {[num2str(f_cut(i)),' Hz']}, 'LineWidth', 2);
    set(get(get(x,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

%%% Comments:
% The SPL vs Octave plot matches the FFT plot. We relate the Weighted plots
% to the non-weighted ones and see that the relative slope of the weighted
% signals to the non-weighted one are preserved. We see that the C-Weighted
% and the non-weighted have very similar slopes to their FFT, so and they
% also almost overlap on the SPL plot. In fact, once the FFT of the
% C-Weighted signal starts to drop off it also drops off in the SPL plot.
% We also consider the A - Weighted signal, where we again see that the SPL
% follows the slope of the FFT in comparison to the non-weighted FFT. It
% begins with a lower slope than the non-weighted, and then the slope
% gets closer and closer to the non-weighted, until it surpasses it. This
% is exactly matched in the SPL plot. 
%% New Functions
%Function inspired by M.Sc. Eng. Hristo Zhivomirov  
function [p] = pinknoise(N,fs)
%[p,P] = pinknoise(N,fs)
% generates N - element pinknoise with sampling frequency fs

X = ones(1,N);
%generate phase
%since we want to decrease power by 1/(sqrt(f)) we can just do that to the
%phase
f = 1:floor(N/2);
phase = exp(-1j*2*pi*randn(1,floor(N/2)))./sqrt(f);
%Set X based on the end (as my ifft is [fs/2 ... 0... ], so if odd we don't
%have fs/2 and if even we do
X(end-(floor(N/2))+1:end) = phase;
X(end - 2*(floor(N/2))+1: end-(floor(N/2))) = conj(phase(end:-1:1));

p = (my_ifft (X,N, fs));
end

%% Functions
%%% my_fft
function [X,f,sum_check_t,sum_check_f] = my_fft(x,N, fs)
%[X,f,sum_check_t,sum_check_f,checkshift] = my_fft(x,N,fs)
%Takes the N point FFT of x and multiplies by 1/fs. Shifts the reference
%output and provides the shifted frequencies for plotting. Two outputs to
%check if the FT was performed properly, and one more to check if the shift
%is like fftshift.
%Make input row vector
x= x(:).';
dt = 1/fs;
df = fs/N;
%Take FT
Xpshift = fft(x,N)*dt;
%Checking if FT was good using Parseval's
sum_check_t= sum(x.^2)*dt;
sum_check_f= (sum(abs(Xpshift).^2))*df;
%shifting X
X = my_fftshift(Xpshift, N);
%Since my shift is like fftshift, the fs/2 component goes to -fs/2, which
%should be the first element of f:
f = [-floor(N/2):-1,0:ceil(N/2)-1]*df;
end
%%% my_ifft
function [x, sum_check_t,sum_check_f] = my_ifft(X,N,fs)
%[x,sum_check_t,sum_check_f] = my_ifft(X,N,fs)
%Takes the N point iFFT of X after dividing X by dt. Two outputs to
%check if the FT was performed properly.
dt = 1/fs;
df = fs/N;
%Unshift - Note that this is NOT my_fftshift
X = [X(floor(N/2)+1:end), X(1:floor(N/2))];
%IFT
x = ifft(X)/dt;
%Checking if FT was good using Parseval's
sum_check_t= sum(x.^2)*dt;
sum_check_f= (sum(abs(X).^2))*df;
end
%%% my_fftshift
function X = my_fftshift(X, N)
%X = myfftshift(X,N)
% My interpertation of fftshift. Takes the input and divides it into 2,
% when the size is odd, we let the make the first half (indices 1 to
% ceil(half)) be the bigger "half", and swap the two sections.
X = [X(ceil(N/2)+1:end), X(1:ceil(N/2))];
end
%%% White noise
function [w] = whitenoise(N,fs)
%[w] = whitenoise(N,fs)
%Generates white noise of size N. 

%Initialize to all 1s
X = ones(1,N);
%generate phase
phase = exp(1j*2*pi*randn(1,floor(N/2)));
%Set X based on the end (as my ifft is [fs/2 ... 0... ], so if odd we don't
%have fs/2 and if even we do
X(end-(floor(N/2))+1:end) = phase;
X(end - 2*(floor(N/2))+1: end-(floor(N/2))) = conj(phase(end:-1:1));

w = (my_ifft (X,N, fs));
end
%%% X2Gxx

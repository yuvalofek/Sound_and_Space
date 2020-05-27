%%% Yuval Epstain Ofek 
%%% Sound and Space - Problem Set 4

%%
clear all; close all; clc;

Recordsize = 500;
PercentOverlap = 75;
%a)
[HW4Spect,fs] = audioread('HW4_spect.wav');
T = size(HW4Spect,1)/fs;

%b)
[Gxx, fg, rec] = plotSpect(HW4Spect, Recordsize, PercentOverlap, fs);
title('Spectrogram of HW4\_spect.wav','FontSize', 16, 'FontWeight', 'bold')
%%% c) Comments
% The spectogram depicts a decrease in frequency over a small time range
% coupled with a sudden pulse in sound intensity followed by a mellowing
% down and return to the previous intensity level and a constant frequency.
% It feels like this is a recording of the doppler effect.We have a number 
% of dominant frequencies at the beginning of the recording, we then see an
% increase in intensity of the peaks as well as a curving down of the lines
% that represented our dominant frequencies, all of which is followed by a 
% decrease in intensity and straightening of the lines at lower frequencies
% compared to before. We can even say that the source was coming towards 
% the microphone and passed the microphone at a bit under 4 seconds. The
% curve represents changing dominant frequencies, in this case decreasing.

%d)
c = 343; %speed of sound in air [m/s]

% The reading is a bit noisy so I ended up getting the values from the
% plot. I looked for the bottom-most horizontal line that had a decent
% intensity. 
f_high = 420;
f_low = 270;

% f1*(1-v/c) = f2*(1+v/c) => v = (f1 - f2)c/(f1+f2)
 v = (f_high - f_low)*c/(f_high  + f_low)

%%% Bonus:
% We end up only getting the radial frequency shift from the microphone,
% which will not really help if the source is moving tangentially or really
% far away from the source:
%       ^
%       |
%       |
%       |  <-source path                                        * <- Mic
%       |
%       |
%       |<----------------------------------------------------->
%       |   Large distance between source path and microphone
%       |
%    
% Doppler shift only accounts for velocity moving in the radial direction. 
% This means that we don't know how the source is moving tangentially to
% the microphone (say the source was moving in circles around the
% microphone, the frequency of the sound it generates at the microphone
% will NOT change). We can account for this by adding a second microphone
% somewhere not by the first (this all should be done relative to the
% source of course) and maybe a third microphone if we want to consider
% vertical motion as well. This way, the speed of the source that is not
% accounted by one microphone (when the source moves tangentially compared
% to that microphone) will be accounted for by a different microphone. 

%% Problem 2
clear all;close all;clc;
%%% a)
[pulsenoise,fs] = audioread('HW2_pulsenoise.wav');

%%% b)
% Auto-correlation
[Rxx, Tau] = xy2Rxy(pulsenoise,pulsenoise,fs);

% Plotting
figure
plot(Tau, Rxx, 'LineWidth', 2)
title('Autocorrelation of HW2\_pulsenoise.wav','FontSize', 16, 'FontWeight', 'bold')
ylabel('Magnitude (WU)','FontSize', 16, 'FontWeight', 'bold')
xlabel('Time (s)','FontSize', 16, 'FontWeight', 'bold')
xlim([min(Tau), max(Tau)])
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Rxx_{pulsenoise}', 'FontSize', 16)
%%% Observations:
% Looking at the plot, we see that we have peaks every 0.08533 (roughly),
% so we can say that our pulse is that long. This corresponds to 1024 
% indices. We also see that there are no gaps in pulsenoise, that the
% signal repeats without breaks for the entirety of pulsenoise. This allows
% us to say that there are around 262196/1024 ~ 256 averages.


%%% c)
pulse = audioread('HW2_pulse.wav');
%Assuming same sampling frequency

%%% d)
[Rxy, Tau2] = xy2Rxy(pulse,pulsenoise,fs);
figure
plot(Tau2, Rxy, 'LineWidth', 2)
title('Cross-Correlation between HW2\_pulse.wav and HW2\_pulsenoise.wav','FontSize', 16, 'FontWeight', 'bold')
ylabel('Magnitude (WU)','FontSize', 16, 'FontWeight', 'bold')
xlabel('Time (s)','FontSize', 16, 'FontWeight', 'bold')
xlim([min(Tau2), max(Tau2)])
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Rxy_{pulsenoise&pulse}', 'FontSize', 16)

%%% Using this...
% We see that the largest maxima is at t = 0.00375 s, so this should be the
% time it takes pulse to show up (as closely as possible) in pulsenoise.
% Now we can take that time and the speed the sound travels in and
% calculate the distance between the source and the microphone:

v_sound = 343; %Speed of sound in air
[~, Imx] = max(Rxy);
tau0 = Tau2(Imx);

d = (v_sound * tau0)


%% New functions
%%% plotSpect
function [Gxx, fg,records] = plotSpect(x, NperChunk, PercentOver, fs)
%%% [Gxx, fg] = plotSpect(x, NperChunk, PercentOver, fs)
%%% Makes a spectogram of the data in x by taking records of length
%%% NperChunk, that overlap PercentOver percent, for a size of x, N, and fs
%%% sampling frequency
N = max(size(x));

T = N/fs;

Nover = round(NperChunk*PercentOver/100);
Nskip = NperChunk-Nover;
records = floor((N-Nover)/Nskip);

%initialize
xch = zeros(records, NperChunk);
Gxx = zeros(records,ceil(NperChunk/2)+1);
fg = zeros(1,ceil(NperChunk/2)+1);

%generate the Gxxs
for i = 1:records
    xch(i,:) = x(1+(i-1)*(Nskip): NperChunk+(i-1)*(Nskip));
    [Gxx(i,:),fg] = x2Gxx(xch(i,:),fs);
end
%make relative to max Gxx
mxGxx = max(max(Gxx));
GxxdB = 20*log(Gxx/mxGxx);

%spectogram
figure
imagesc('xdata',(1:records)*T/records, 'ydata',fg,'Cdata',GxxdB.')
xlim([0,T])
ylim([0,fs/2])
c = colorbar('Fontsize', 16);
c.Label.String = 'Magnitude (dB re: max)';
xlabel('Time (s)','FontSize', 16, 'FontWeight', 'bold')
ylabel('Frequency (Hz)','FontSize', 16, 'FontWeight', 'bold')
colormap jet
end
%%% xy2CSD
function [Gxy, fc] = xy2CSD (x,y, fs)
N = max(size(x));
T = N/fs;
[X,f] = my_fft(x,N,fs);
[Y,~] = my_fft(y,N,fs);

[Gxy,fc] = XY2CSD(X,Y,f,T,fs);

end
%%% XY2CSD
function [Gxy, fc] = XY2CSD (X,Y,f,T,fs)
Sxy = (Y.*X'.')/T;
I0 = find(f == 0);
Gxy = [Sxy(I0), 2*Sxy(I0+1:end), Sxy(1)];
fc = [f(I0:end), fs/2];
end
%%% XY2Rxy
function [Rxy, Tau] = XY2Rxy(X,Y,T,fs)
Sxy = (Y.*conj(X))/T;
N = max(size(X));
fullRxy = my_ifft(Sxy,N,fs);
Rxy = my_fftshift(fullRxy, N);
Tau = [-floor(N/2):-1,0:ceil(N/2)-1]/fs;
end
%%% xy2Rxy
function [Rxy, Tau] = xy2Rxy(x,y,fs)
N = max(size(x));
T = N/fs;
[X,~] = my_fft(x,N,fs);
[Y,~] = my_fft(y,N,fs);
[Rxy, Tau] = XY2Rxy(X,Y,T,fs);
end


%% Functions:
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
%%% X2Gxx
function [Gxx,fg] = X2Gxx(X,f,T,fs)
%[Gxx,fg] = getGxx(X,f,T,fs)
%Given a linear spectrum, T for the sample, and the sampling frequency,
%calculates the Gxx of the linear spectrum
Sxx = (X.*conj(X))/T;
I0 = find(f == 0);
Gxx = [Sxx(I0), 2*Sxx(I0+1:end), Sxx(1)];
fg = [f(I0:end), fs/2];
end
%%% x2Gxx
function [Gxx,fg] = x2Gxx(x,fs)
% [Gxx,fg] = x2Gxx(x,fs)
% From time series calculates Gxx. 
N = max(size(x));
T = N/fs;
[X,f] = my_fft(x,N,fs);
[Gxx,fg] = X2Gxx(X,f,T,fs);
end


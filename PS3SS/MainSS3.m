% Yuval Epstain Ofek
% Sound & Space HW3
%%
clear all;close all;clc;
%%% Problem 1: Created functions PSD_Asynch & PSD_Synch at the bottom of
%%% the document

%%% a)
[pulsenoise,fs] = audioread('HW2_pulsenoise.wav');
N = size(pulsenoise, 1);
T = N/fs;
tn = (0:(N-1))/fs;

%%% b)
[Gxx,fg] = x2Gxx(pulsenoise,fs);
mxG = max(Gxx);
Gxx_dB = 10*log(Gxx/mxG);

figure
plot (fg,Gxx_dB, 'LineWidth' , 1.5)
title('Gxx of pulsenoise.wav', 'FontSize', 18, 'FontWeight', 'bold')
xlabel('Frequency(Hz)','FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude (dB re: max)','FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx', 'FontSize', 16)
xlim([0,fs/2])
ylim([min(Gxx_dB)-5, max(Gxx_dB)+5])

%%% Comments:
% It is hard to determine the number of overtones. We see at least 7 peaks
% above the noise, and at most 11 (but some of these are really close to
% the level of the noise, so it is hard to know if they are truly overtones
% or actually just noise). If we talk about reliability, then the number of
% overtones we can clearly distinguish is 7.

%%% c)
[pulse] = audioread('HW2_pulse.wav');        %since fs is the same we ignore it here
Np = size(pulse, 1);
Tp = N/fs;
t = (0:(Np-1))/fs;

figure
plot(t,pulse)
title('time series of pulse.wav', 'FontSize', 18, 'FontWeight', 'bold')
xlabel('time (s)','FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude (WU)','FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Time series', 'FontSize', 16)
xlim([0,t(end)])
ylim([min(pulse)-0.5, max(pulse)+0.5])

%We look at the time series of pulse and see that it is zero until index
%200, proceeds with a square pulse for 400 indices, and further has 424
%indices at zero following the pulse. This allows us to see that Nsamp
%should be 400, Ngap should be 624 (424+200), and the Nstart should be 200.
%We consider this problem, and assume that what this scenario models is us
%starting to record and sending a series of pulses. This means that there
%should be some delay to get the noise that we sent, so we increase nstart
%by a bit to 210. 

[Gxxout,fgo] = PSD_Synch(pulsenoise,fs, 2, 400 ,210, 624);
mxGo = max(Gxxout);
Gxxo_dB = 10*log(Gxxout/mxGo);

figure
plot(fgo,Gxxo_dB, 'LineWidth', 2)
title('Gxx of pulsenoise.wav - Time AVG', 'FontSize', 18, 'FontWeight', 'bold')
xlabel('Frequency(Hz)','FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude (dB re: max)','FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx_{timeAVG}', 'FontSize', 16)
xlim([0,fs/2])
ylim([min(Gxxo_dB)-5, max(Gxxo_dB)+5])

%%% Comments:
% Now we can clearly see 11 overtones. This is considerably more than we
% previously saw with just RMS.

%%% d)
figure 
plot( fg,Gxx_dB,fgo,Gxxo_dB, 'LineWidth', 2)
title('Comparing RMS and Time Averaging', 'FontSize', 18, 'FontWeight', 'bold')
xlabel('Frequency(Hz)','FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude (dB re: max)','FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx_{TimeAVG}', 'Gxx_{RMS}', 'FontSize', 16)
xlim([0,fs/2])
ylim([min(Gxx_dB)-10, max(Gxx_dB)+10])

%%% Comments:
% We notice far more overtones in the synchronous time average compared to
% the RMS. We also see that the magnitude for these overtones (in dB re:max)
% is slightly greater in the time average than in the RMS. If we consider 
% the rest of the Gxx, we see that the RMS is more jagged than the time
% average (goes up and down more frequently and by greater magnitudes).
% This shows us that the noise variance was reduced going from RMS to time
% averaging. If we draw a line along the top of the two Gxx plots, we see
% that in the two effects together allowed for more significant peaks
% and troughs in the time averaging Gxx, and so more overtones can be
% noticed. We also see that the "data", or memory, that the time average 
% takes to store is less than the RMS(comparing the two, the RMS Gxx 
% requires 131099 elements while the time average only needs 201, this is 
% a considerable difference!), so not only do we get more meaningful
% results, we can store it easier.

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
Sxx = (X.*X'.')/T;
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
%%% PSD_Asynch- asynchronous PSD analyses 
function [Gxxout,fg] = PSD_Asynch( x, fs, Num_avg, PSD_type)
% [Gxxout,fg] = get_PSD_Asynch( x, fs, Num_avg, PSD_type)
% PSD_type = 1 - RMS
%          = 2 - Linear Averaging
%          = 3 - time averaging
% Calculates the PSD of the time sample in the way described by PSD_type.
% For PSD_type = 1, takes the full time series PSD, for PSD_type = 2, takes
% an average Gxx of the ffts taken after cropping the time series into
% Num_avg parts. For PSD_type = 3, takes the time series, crops it into
% Num_avg parts, takes their average, and then computes the Gxx of that.

if PSD_type == 1   %Pure RMS
    [Gxxout, fg] = x2Gxx(x,fs);

elseif PSD_type == 2    %Linear averaging
    N = max(size(x));
    %number of samples per chunk
    Nperch = floor(N/Num_avg);
    %initialize
    xch = zeros(Num_avg, Nperch);
    Gxx = zeros(Num_avg,ceil(Nperch/2)+1);
    %generate the Gxxs
    for i = 1:Num_avg
        %split the time series into chunks, take fft of each, and compute
        %Gxx
        xch(i,:) = x((i-1)*Nperch+1:Nperch+(i-1)*Nperch);
        [Gxx(i,:), fg(:)] = x2Gxx(xch(i,:),fs);
    end
    %output the mean Gxx
    Gxxout = mean(Gxx,1);
    
else %time averaging
    N = max(size(x));
    %number of samples per chunk
    Nperch = floor(N/Num_avg);
    %initialize
    xch = zeros(Num_avg, Nperch);
    %split time series
    for i = 1:Num_avg
        xch(i,:) = x((i-1)*Nperch+1:Nperch+(i-1)*Nperch);
    end
    %average, take fft, and fing Gxx
    xnew = mean(xch,1);
    [Gxxout,fg] = x2Gxx(xnew,fs); 
end
end
%%% PSD_Synch - with index parameters to allow for synchronous averaging
function [Gxxout,fg] = PSD_Synch( x,fs, PSD_type, Nsamp,Nstart, Ngap)
% [Gxxout,fg] = get_PSD_Synch( x, fs, PSD_type, timevector)
% PSD_type = 1 - linear averaging
%          = 2 - time averaging
%Takes the averaged PSD of x with a given sampling frequency. To pick the
%time series chunks to be averaged, we set parameters: Nsamp,Nstart,Ngap,and
%periods. Nsamp is the number of samples per time series chunk, Nstart is
%the offset from the start of the time series to start taking our time
%series chunk, Ngap is the time break between the time series chunks, and
%periods is the number of time series chunks that are collected.

%Splitting x
N = max(size(x));
Nskip = Ngap+Nsamp;
periods = floor(N/Nskip);
xch = zeros(periods, Nsamp);
for i = 1:periods
    xch(i,:) = x(Nstart+(i-1)*(Nskip): Nstart+Nsamp-1+(i-1)*(Nskip)).';
end
if PSD_type == 1
    %initialize
    Gxx = zeros(periods,ceil(Nsamp/2)+1);
    %generate the Gxxs
    for i = 1:periods
        [Gxx(i,:),fg(:)] = x2Gxx(xch(i,:),fs);
    end
    Gxxout = mean(Gxx,1);
else 
    xnew = mean(xch,1);
    [Gxxout, fg] = x2Gxx(xnew,fs);    
end
end
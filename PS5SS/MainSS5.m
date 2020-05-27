%Yuval Epstain Ofek
%S&S5

%% Problem 1
clear all;close all;clc;

Plength = 50000;
Rlength = Plength+1000;
numAvg = 2;
fWave = 2;
fs = 12000;
f1 = 1000;
f2 = 5000;
Ntot = Rlength*numAvg;

%%% I change some of the frequencies when calling the code to make it
%%% easier to interprate the graphs.
% Impulse
[y] = GenExc(Rlength,numAvg);
figure 
subplot(3,2,1)
plot((0:Ntot-1)/fs, y, 'LineWidth', 1.5)
title('Impulse','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Impulse', 'FontSize', 16)
xlim([0,(Ntot-1)/fs])
ylim([-0.1, 1.1])

% White noise
[y] = GenExc(Rlength,numAvg,Plength,fs);
subplot(2,2,2)
plot((0:Ntot-1)/fs, real(y), 'LineWidth', 1)
title('White Noise','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('White Noise', 'FontSize', 16)
xlim([0,(Ntot-1)/fs])
ylim([-250, 250])

% CW
[y] = GenExc(Rlength,numAvg,Plength,fs,f1/100);
subplot(3,2,3)
plot((0:Ntot-1)/fs, y, 'LineWidth', 1)
title('CW Pulse','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('CW Pulse', 'FontSize', 16)
xlim([0,(Ntot-1)/fs])
ylim([-1.1, 1.1])

% Linear Chirp - reduced f1 & f2 so it is easy to see
[y] = GenExc(Rlength,numAvg,Plength,fs,f1/100,f2/100,1);
subplot(2,2,4)
plot((0:Ntot-1)/fs, y, 'LineWidth', 1)
title('Linear Chirp','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Linear Chirp', 'FontSize', 16)
ylim([-1.1, 1.1])

% Logarithmic Chirp
[y] = GenExc(Rlength,numAvg,Plength,fs,f1*10,f2*100,0);
subplot(3,2,5)
plot((0:Ntot-1)/fs, y, 'LineWidth', 1.5)
title('Logarithmic Chirp','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Logarithmic Chirp', 'FontSize', 16)
xlim([0,(Ntot-1)/fs])
ylim([-1.1, 1.1])


%%% Comments:
% The excitations look as expected. We see that the excitation is repeated
% exactly on the second recorda and that the record length can be more than
% the pulselength. 
%% Problem 2
clear all;close all;clc;
T = 10;
Tp = 5;
fs = 12000;
f1 = 2000;
f2 = 5000;

Plength = Tp*fs;
Rlength = T*fs;
numAvg = 5;

ExciteTest(Rlength,numAvg,Plength,fs,f1,f2)
Tp = 10;
Plength = Tp*fs;
ExciteTest(Rlength,numAvg,Plength,fs,f1,f2)


%%% Comments (a)
% Rxx - The Rxx for the impulse and white noise did not change at all. This
% makes sense as the number of records didn't change and the two
% excitations are distinct across the entire pulse length. The shape of the
% Rxx of the linear chirp didn't change, but the magnitude increased when
% Tp was increased. This makes sense, as now the excitation changes more
% slowly, so the Rxx should increase. The shape of the CW pulse changed
% from a series of diamonds to that of a rectangular block. This is
% reasonable, as before the pulse took the time for half a record and after
% it was the entire record. As for the logarithmic chirp, the Rxx slightly
% changed in shape with the increase in Tp, making the side 'lobes' more
% pronounced and distinct. The amplitude also increased. 
% Gxx - The Gxx of the impulse didn't change. The Gxx for the sinusoidal
% and the linear chirp seemed to look the same but with quadruple and
% double the magnitude respectively. The Gxx of the logarithmic pulse also
% looked the same with an increased magnitude. The white noise Gxx for the
% smaller Tp looked jagged and similar to white noise itself, but after
% doubling the Tp, the Gxx became a magnitude 1 flat spectrum (plus little
% noise). 


%%% Comments (b)
% One type of excitation might be more preferable over another based on the
% test we want to run. If we want to test 1 frequency, we better us a
% sinusoidal excitation. For a range, we can try to discretize the range
% using CW pulses, or use one of the chirps. If we just was to find how the
% system reacts to all frequencies, we should use white noise or impulses
% (but impulses tend to not work best in real life scenarios). Given what
% we want to find out about our system, different excitations will work
% better.

%% Problem 3
clear all;close all; clc;
% a)
T =1;
Tp = 0.1;
fs = 12000;
f1 = 2000;

Plength = Tp*fs;
Rlength = T*fs;
numAvg = 1;

t = (0:Rlength-1)/fs;
[y] = GenExc(Rlength,numAvg,Plength,fs,f1);
[Gyy,fy] = x2Gxx(y,fs);

%b)
Lhan = Plength/2;
han = hann(Lhan);
hanwin = [han(1:Lhan/2);ones(Plength-Lhan-1,1);han(Lhan/2:end)]; 
ywin = zeros(Rlength,1);
ywin(1:Plength) = hanwin.*y(1:Plength);

plot(t,y,t,ywin,'LineWidth', 1.5)
title('Time domain of y','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('y','windowed y', 'FontSize', 16)


[Gywyw,fyw] = x2Gxx(ywin,fs);

figure
subplot(2,1,1)
plot(fy,Gyy,'LineWidth', 1.5)
title('Original Gyy','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gyy', 'FontSize', 16)

subplot(2,1,2)
plot(fyw,Gywyw,'LineWidth', 1.5)
title('Gyy of Windowed y','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gyy', 'FontSize', 16)

%%% Comments
% The taper makes the lowers the magnitude of the PSD. On the other hand it
% also reduces the number of side-lobes that we see in the Gxx plot.



%% Problem 4
clear all;close all; clc;

% a)
%Parameters
T =1;
Tp = 0.1;
fs = 12000;
f1 = 100;
f2 = 5000;

Plength = Tp*fs;
Rlength = T*fs;
numAvg = 1;

%generate output
t = (0:Rlength-1)/fs;
[ylin] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,1);
[ylog] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,0);

%windowing
Lhan = Plength/2;
han = hann(Lhan);
hanwin = [han(1:Lhan/2);ones(Plength-Lhan-1,1);han(Lhan/2:end)]; 

ylinwin = zeros(Rlength,1);
ylogwin = zeros(Rlength,1);

ylinwin(1:Plength) = hanwin.*ylin(1:Plength);
ylogwin(1:Plength) = hanwin.*ylog(1:Plength);

%plot results
figure 
subplot(2,1,1)
plot(t,ylin,t, ylinwin,'LineWidth', 1.5)
title('Linear Chirp','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('y','windowed y', 'FontSize', 16)

subplot(2,1,2)
plot(t,ylog,t, ylogwin,'LineWidth', 1.5)
title('Logarithmic Chirp','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('y','windowed y', 'FontSize', 16)

%Gxxs
[Gyylin,fylin] = x2Gxx(ylin,fs);
[Gywywlin,fywlin] = x2Gxx(ylinwin,fs);
[Gyylog,fylog] = x2Gxx(ylog,fs);
[Gywywlog,fywlog] = x2Gxx(ylogwin,fs);


figure
subplot(2,1,1)
plot(fylin,Gyylin,'LineWidth', 1.5)
hold on
plot(fywlin,Gywywlin,'LineWidth', 1.5)

title('Original Gyy','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Original Gyy','Windowed Gyy', 'FontSize', 16)

subplot(2,1,2)
plot(fylog,Gyylog,'LineWidth', 1.5)
hold on
plot(fywlog,Gywywlog,'LineWidth', 1.5)
title('Gyy of Windowed y','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gyy','Windowed Gyy', 'FontSize', 16)

%%% Comments:
% The PSD plots became a lot smoother than before, but we lost the edges of
% our the PSD peak and also lost some of the magnitude.

%% d)
%repeat above for different f2
T =1;
Tp = 0.1;
fs = 12000;
f1 = 100;
f2 = 8000;

Plength = Tp*fs;
Rlength = T*fs;
numAvg = 1;

t = (0:Rlength-1)/fs;
[ylin] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,1);
[ylog] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,0);

Lhan = Plength/2;
han = hann(Lhan);
hanwin = [han(1:Lhan/2);ones(Plength-Lhan-1,1);han(Lhan/2:end)]; 

ylinwin = zeros(Rlength,1);
ylogwin = zeros(Rlength,1);

ylinwin(1:Plength) = hanwin.*ylin(1:Plength);
ylogwin(1:Plength) = hanwin.*ylog(1:Plength);

figure 
subplot(2,1,1)
plot(t,ylin,t, ylinwin,'LineWidth', 1.5)
title('Linear Chirp','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('y','windowed y', 'FontSize', 16)

subplot(2,1,2)
plot(t,ylog,t, ylogwin,'LineWidth', 1.5)
title('Logarithmic Chirp','FontSize', 16)
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('y','windowed y', 'FontSize', 16)


[Gyylin,fylin] = x2Gxx(ylin,fs);
[Gywywlin,fywlin] = x2Gxx(ylinwin,fs);
[Gyylog,fylog] = x2Gxx(ylog,fs);
[Gywywlog,fywlog] = x2Gxx(ylogwin,fs);


figure
subplot(2,1,1)
plot(fylin,Gyylin,'LineWidth', 1.5)
hold on
plot(fywlin,Gywywlin,'LineWidth', 1.5)

title('Original Gyy','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Original Gyy','Windowed Gyy', 'FontSize', 16)

subplot(2,1,2)
plot(fylog,Gyylog,'LineWidth', 1.5)
hold on
plot(fywlog,Gywywlog,'LineWidth', 1.5)
title('Gyy of Windowed y','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude [WU]', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gyy','Windowed Gyy', 'FontSize', 16)

%%% Comments:
% When we let f2 be 8000 > fs/2, we get some really funky results. We
% significant power at the higher frequencies (close to fs/2) in the
% linear case oscillates and grows for both the windowed and not, and for
% the logarithmic case the power does increase at the higher end, but not
% as much as the linear chirp, and this is not represented in the windowed
% PSD plot. 
%% New Functions
%%% White noise
function [w] = whitenoise(N,fs)
%[w] = whitenoise(N,fs)
%Generates white noise of size N. 

%Initialize to all 1s
X = ones(1,N);
%generate phase
phase = exp(-1j*2*pi*randn(1,floor(N/2)));
%Set X based on the end (as my ifft is [fs/2 ... 0... ], so if odd we don't
%have fs/2 and if even we do
X(end-(floor(N/2))+1:end) = phase;
X(end - 2*(floor(N/2))+1: end-(floor(N/2))) = conj(phase(end:-1:1));

w = real(my_ifft (X,N, fs));
end
%%% Excitation generator
function [x] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,isLinear)
%[x] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,isLinera)
% Excitations: 1 - impulse
%              2 - CW Wave (sinusoidal)
%              3 - Linear Chirp
%              4 - Logarithmic Chirp
%              5 - White Noise
% [x] = GenExc(Rlength,numAvg) 
%       Returns a vector with numAvg pulses spaced by Rlength
% [x] = GenExc(Rlength,numAvg,Plength, fs)
%       Returns a vector with white-noise of length Plength spaced by
%       Rlength. The white noise is the SAME each time. 
% [x] = GenExc(Rlength,numAvg,Plength,f1)
%       Returns a vector with sinosoidal of frequency f1 of length Plength
%       and spaced by Rlength
% [x] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,isLinear)
%       Returns a vector with chirp signal going from f1 to f2 of length
%       Plength and spaced by Rlength. Chirp is linear id isLinear is 1 and
%       logarithmic id isLinear is anything else.

%select type based on number of arguments
if nargin == 2
    extype = 1; %only Rlength and numAvg
elseif nargin == 4
    extype = 5; 
elseif nargin == 5
    extype = 2;
else
    if (isLinear  == 1)
        extype = 3;
    else 
        extype = 4;
    end
end    
        
%total length of output vector
Tlength = Rlength*numAvg;
x = zeros(Tlength,1);
switch extype
    case 1
        %impulse
        for ii = 1:numAvg
            x(Rlength*(ii-1)+1) = 1;
        end
    case 2
        %sinusoidal
        t = (0:Plength-1)/fs;
        sn = sin(2*pi*f1*t);
        for ii = 1:numAvg
            x(Rlength*(ii-1)+1:Plength+Rlength*(ii-1)) = sn;
        end
    case 3
        Tp = Plength/fs;
        %linear chirp
        t = (0:Plength-1)/fs;
        phi = pi*(f2-f1)/Tp*t.^2 + 2*pi*f1*t;
        s = sin(phi);
        for ii = 1:numAvg
            x(Rlength*(ii-1)+1:Plength+Rlength*(ii-1)) = s;
        end
    case 4
        %logarithmic chirp
        Tp = Plength/fs;
        t = (0:Plength-1)/fs;
        phi = 2*pi*f1*Tp/(log(f2/f1))*(f2/f1).^(t/Tp)-2*pi*f1*Tp/(log(f2/f1));
        s = sin(phi);
        for ii = 1:numAvg
            x(Rlength*(ii-1)+1:Plength+Rlength*(ii-1)) = s;
        end
    case 5
        wn = whitenoise (Plength,fs);
        for ii = 1:numAvg
            x(Rlength*(ii-1)+1:Plength+Rlength*(ii-1)) = wn;
        end
end
end
%%% Test Excitations
function ExciteTest(Rlength,numAvg,Plength,fs,f1,f2)

%Pulse
[y] = GenExc(Rlength,numAvg);
[Gyyp,fgp] = x2Gxx(y,fs);
[Ryyp, Taup] = xy2Rxy(y,y,fs);

% White noise
[y] = GenExc(Rlength,numAvg,Plength,fs);
[Gyyw,fgw] = x2Gxx(y,fs);
[Ryyw, Tauw] = xy2Rxy(y,y,fs);
%CW
[y] = GenExc(Rlength,numAvg,Plength,fs,f1);
[GyyCW,fgCW] = x2Gxx(y,fs);
[RyyCW, TauCW] = xy2Rxy(y,y,fs);

%Chirp Linear
[y] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,1);
[Gyyclin,fgclin] = x2Gxx(y,fs);
[Ryyclin, Tauclin] = xy2Rxy(y,y,fs);

%Chirp Log
[y] = GenExc(Rlength,numAvg,Plength,fs,f1,f2,0);
[Gyyclog,fgclog] = x2Gxx(y,fs);
[Ryyclog, Tauclog] = xy2Rxy(y,y,fs);

figure 
sgtitle('Gxx','FontSize', 20)

subplot(3,2,1)
plot(fgp, Gyyp, 'LineWidth', 1.5)
title('Pulse','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx', 'FontSize', 16)
ylim([1.1*min(Gyyp), 1.1*max(Gyyp)])

subplot(2,2,2)
plot(fgw, Gyyw, 'LineWidth', 1.5)
title('White Noise','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx', 'FontSize', 16)
ylim([1.1*min(Gyyw), 1.1*max(Gyyw)])


subplot(3,2,3)
plot(fgCW,GyyCW, 'LineWidth', 1.5)
title('Sinusoidal','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx', 'FontSize', 16)
ylim([1.1*min(GyyCW), 1.1*max(GyyCW)])


subplot(2,2,4)
plot(fgclin,Gyyclin, 'LineWidth', 1.5)
title('Linear Chirp','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx', 'FontSize', 16)
ylim([1.1*min(Gyyclin), 1.1*max(Gyyclin)])

subplot(3,2,5)
plot(fgclog,Gyyclog, 'LineWidth', 1.5)
title('Logarithmic Chirp','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Gxx', 'FontSize', 16)
ylim([1.1*min(Gyyclog), 1.1*max(Gyyclog)])

%%%%%
figure 
sgtitle('Rxx','FontSize', 20)

subplot(3,2,1)
plot(Taup, Ryyp, 'LineWidth', 1.5)
title('Pulse','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Rxx', 'FontSize', 16)
ylim([1.1*min(Ryyp), 1.1*max(Ryyp)])

subplot(2,2,2)
plot(Tauw, Ryyw, 'LineWidth', 1.5)
title('White Noise','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Rxx', 'FontSize', 16)

subplot(3,2,3)
plot(TauCW,RyyCW, 'LineWidth', 1.5)
title('Sinusoidal','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Rxx', 'FontSize', 16)

subplot(2,2,4)
plot(Tauclin,Ryyclin, 'LineWidth', 1.5)
title('Linear Chirp','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Rxx', 'FontSize', 16)

subplot(3,2,5)
plot(Tauclog,Ryyclog, 'LineWidth', 1.5)
title('Logarithmic Chirp','FontSize', 16)
xlabel('frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Amplitude', 'FontSize', 16, 'FontWeight', 'bold')
grid on;
grid minor
ax = gca;
ax.GridAlpha = 0.5;
ax.FontSize = 16;
legend('Rxx', 'FontSize', 16)
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
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


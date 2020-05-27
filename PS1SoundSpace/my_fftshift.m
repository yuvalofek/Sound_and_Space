function X = my_fftshift(X, N)
%X = myfftshift(X,N)
% My interpertation of fftshift. Takes the positive and zero indices
%(assuming that zero is centered) and switches with the negative ones. 
X = [X(ceil(N/2)+1:end), X(1:ceil(N/2))];
end
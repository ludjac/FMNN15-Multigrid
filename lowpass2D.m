function rf = lowpass2D(rf)
% lowpass filter which dampends high frequencies

rf = conv2(rf,[1 2 1; 2 4 2; 1 2 1]/16,'same');

end
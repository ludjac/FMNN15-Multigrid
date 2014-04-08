function ef = FMGprolong2D(ec)
% prolongs from coarse to fine grid
b = zeros(size(ec)*2 + 1);
b(2:2:end,2:2:end) = ec;
ef = conv2(b, [1 2 1; 2 4 2; 1 2 1]/4, 'same');
end
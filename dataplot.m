figure(i)
for i=1:1:6
    subplot(3, 2, i);
    ufinal = plotdata{i};
    h = surf(X,Y,ufinal, 'facecolor','interp', 'edgecolor','none');
    
    %Controls the colors of the drawn mesh
    set(h,'CData',32*ufinal+32);
    
    %Change the height data to be that of the new time step
    set(h,'ZData',ufinal);
end
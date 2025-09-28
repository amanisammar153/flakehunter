function val = umPerPixel(magIndex)
%UMPERPIXEL micrometers per pixel for the given magnification index.
% Matches Utils/conversion_functions.py -> MICROMETER_PER_PIXEL
% Index: 1:2.5x, 2:5x, 3:20x, 4:50x, 5:100x
mp = containers.Map( ...
    {1,     2,      3,      4,      5}, ...
    {3.0754,1.5377, 0.3844, 0.1538, 0.0769} );
assert(isKey(mp, magIndex), 'Unsupported magnification index: %d', magIndex);
val = mp(magIndex);
end

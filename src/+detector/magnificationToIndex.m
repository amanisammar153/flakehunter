function idx = magnificationToIndex(magnification)
%MAGNIFICATIONTOINDEX 2.5->1, 5->2, 20->3, 50->4, 100->5
switch magnification
    case 2.5,  idx = 1;
    case 5,    idx = 2;
    case 20,   idx = 3;
    case 50,   idx = 4;
    case 100,  idx = 5;
    otherwise, error('Unsupported magnification: %g', magnification);
end
end

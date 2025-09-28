function files = listImages(folder)
%LISTIMAGES return full paths of image files in folder (non-recursive).
ex = ["*.png","*.jpg","*.jpeg","*.tif","*.tiff"];
L = [];
for i=1:numel(ex)
    L = [L; dir(fullfile(folder, ex(i)))]; %#ok<AGROW>
end
files = fullfile({L.folder}, {L.name});
end

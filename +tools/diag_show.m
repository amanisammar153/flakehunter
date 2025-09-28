function diag_show(resultsDir, gtCsv, imageBase)
% tools.diag_show('results', 'gt/gt_boxes.csv', '10_x20.png')
% Shows GT (magenta) vs detections (green) on the ORIGINAL image.

assert(isfolder(resultsDir), 'Results folder not found');
assert(isfile(gtCsv), 'GT CSV not found');
assert(ischar(imageBase) || isstring(imageBase), 'imageBase must be like "10_x20.png"');

% --- load GT rows
GT = readtable(gtCsv);
gtRows = GT(strcmpi(string(GT.image), string(imageBase)), :);

% --- find JSON that corresponds to this image (by normalized base)
base = normBase(imageBase);
J = dir(fullfile(resultsDir,'json','*_flakes.json'));
S = []; found = false;
for k=1:numel(J)
    Ss = jsondecode(fileread(fullfile(J(k).folder,J(k).name)));
    key = "";
    if isfield(Ss,'summary')
        if isfield(Ss.summary,'image_rel_path') && ~isempty(Ss.summary.image_rel_path)
            key = string(Ss.summary.image_rel_path);
        elseif isfield(Ss.summary,'image') && ~isempty(Ss.summary.image)
            key = string(Ss.summary.image);
        elseif isfield(Ss.summary,'image_base') && ~isempty(Ss.summary.image_base)
            key = string(Ss.summary.image_base);
        end
    end
    if key == "", key = string(J(k).name); end
    if normBase(key) == base, S = Ss; found = true; break; end
end
assert(found, "No JSON found for " + string(imageBase));

% --- locate the original image
Ipath = '';
if isfield(S.summary,'image') && ~isempty(S.summary.image) && isfile(S.summary.image)
    Ipath = S.summary.image;
elseif isfield(S.summary,'image_rel_path') && ~isempty(S.summary.image_rel_path)
    cand = fullfile(pwd, char(S.summary.image_rel_path));
    if isfile(cand), Ipath = cand; end
end
if isempty(Ipath)
    F = dir(fullfile(pwd,'**', char(imageBase)));
    assert(~isempty(F), 'Could not locate original image: ' + string(imageBase));
    Ipath = fullfile(F(1).folder, F(1).name);
end
I = imread(Ipath);

% --- detections -> [x y w h]
detB = zeros(0,4);
if isfield(S,'flakes') && ~isempty(S.flakes)
    db = vertcat(S.flakes.bbox);
    if ~isempty(db), detB = [[db.x].' [db.y].' [db.w].' [db.h].']; end
end

% --- draw
figure('Name','GT (magenta) vs Detections (green)','NumberTitle','off');
imshow(I); hold on;
for i=1:size(detB,1)
    rectangle('Position', detB(i,:), 'EdgeColor','g', 'LineWidth',1);
end
for i=1:height(gtRows)
    rectangle('Position', [gtRows.x(i) gtRows.y(i) gtRows.w(i) gtRows.h(i)], ...
              'EdgeColor','m', 'LineWidth',2);
end
title(sprintf('%s | det=%d, GT=%d', imageBase, size(detB,1), height(gtRows)));
hold off;
end

function s = normBase(pathOrName)
p = char(string(pathOrName));
[~,nm,~] = fileparts(p); nm = lower(nm);
suffixes = ["_flakes","_overlay","_mask","_gt","_gt_mask"];
changed = true;
while changed
    changed = false;
    for t = suffixes
        if endsWith(nm, t)
            nm = extractBefore(nm, strlength(nm)-strlength(t)+1);
            changed = true;
        end
    end
end
s = string(nm);
end

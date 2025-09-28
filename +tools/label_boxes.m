function label_boxes(imageDir, outCsv)
% LABEL_BOXES  Quick box labeler: draw rectangles on images, save to CSV.
% Usage:
%   tools.label_boxes('flakesSET1', fullfile('gt','gt_boxes.csv'))

if nargin<1 || strlength(imageDir)==0, imageDir = fullfile(pwd,'flakesSET1'); end
if nargin<2 || strlength(outCsv)==0,  outCsv  = fullfile(pwd,'gt','gt_boxes.csv'); end
if ~exist(fileparts(outCsv),'dir'), mkdir(fileparts(outCsv)); end

% List images (robust if +detector helpers are missing)
L = [];
try
    if exist('detector.listImages','file')
        L = detector.listImages(imageDir);
    end
catch, end
if isempty(L)
    exts = ["*.png","*.jpg","*.jpeg","*.tif","*.tiff"];
    for e = exts
        L = [L; string(fullfile(imageDir, e))]; %#ok<AGROW>
    end
    L = arrayfun(@(p) dir(p), L, 'uni',0); L = vertcat(L{:});
    L = arrayfun(@(d) fullfile(d.folder, d.name), L, 'uni',0);
end
if isempty(L), error('No images in %s', imageDir); end

% natural sort if available
try, L = detector.natsortfiles(L); catch, L = sort(L); end

% load existing CSV (append mode)
rows = {};
if isfile(outCsv)
    T0 = readtable(outCsv);
    for i=1:height(T0)
        rows(end+1,:) = {char(T0.image(i)), T0.x(i), T0.y(i), T0.w(i), T0.h(i)}; %#ok<AGROW>
    end
end

hasDrawRect = ~isempty(which('drawrectangle'));
hFig = figure('Name','Label boxes: ←/→, s=save, d=undo, q=quit','NumberTitle','off');
i = 1;
while i >= 1 && i <= numel(L)
    I = imread(L{i});
    clf; imshow(I); title(sprintf('[%d/%d] %s — draw rects; Enter/double-click to close', ...
        i, numel(L), basename(L{i})));

    imgBase = [basename(L{i}) '.png'];
    rects = {}; hs = [];

    keep = true;
    while keep
        if hasDrawRect
    % Make the ROI non-deletable to avoid accidental invalidation
    r = drawrectangle('InteractionsAllowed','all','Deletable',false);
    if isempty(r) || ~isvalid(r), break; end

    pos = [];
    try
        % User double-clicks or presses Enter to finish
        wait(r);
        if isvalid(r)
            pos = round(r.Position);   % [x y w h]
        end
    catch
        % e.g., figure closed or ROI deleted → pos stays empty
    end

    % If ROI was cancelled/deleted, skip cleanly
    if isempty(pos) || any(isnan(pos))
        if isvalid(r), delete(r); end
        % If no boxes yet on this image, offer to skip the image
        if isempty(rects)
            sel = questdlg('No boxes added. Skip this image?', 'Skip', 'Yes','No','Yes');
            if strcmp(sel,'Yes')
                keep = false;   % exit draw loop -> go to navigation menu
                break;
            else
                continue;       % keep drawing
            end
        else
            continue;           % we already have some boxes; just skip this attempt
        end
    end

    hs(end+1) = r; %#ok<AGROW>

else
    % Fallback for older MATLAB without drawrectangle
    try
        h = imrect;
        if isempty(h) || ~isvalid(h), break; end
        p = wait(h);                       % [x y w h]
        if isempty(p) || any(isnan(p))
            if isvalid(h), delete(h); end
            if isempty(rects)
                sel = questdlg('No boxes added. Skip this image?', 'Skip', 'Yes','No','Yes');
                if strcmp(sel,'Yes')
                    keep = false; break;
                else
                    continue;
                end
            else
                continue;
            end
        end
        pos = round(p);
        hs(end+1) = h; %#ok<AGROW>
    catch
        % user cancelled → skip
        continue;
    end
end

rects{end+1} = pos; %#ok<AGROW>


       



        
%%
        sel = questdlg('Add another rectangle?', 'More?', 'Yes','Undo last','No','Yes');
        switch sel
            case 'Undo last'
                if ~isempty(hs), delete(hs(end)); hs(end)=[]; rects(end)=[]; end
            case 'No'
                keep = false;
        end
    end

    % append rectangles for this image
    for k=1:numel(rects)
        r = rects{k};
        rows(end+1,:) = { imgBase, r(1), r(2), r(3), r(4) }; %#ok<AGROW>
    end

    % quick actions (save or undo last saved row for this image)
nav = questdlg('Next action?', 'Navigation','Next','Undo last label','Save','Next');
switch nav
    case 'Undo last label'
        idx = find(strcmp(rows(:,1), imgBase));
        if ~isempty(idx), rows(idx(end),:) = []; end
    case 'Save'
        writeRows(outCsv, rows);
end

% where to move now?
choice = menu('Move to which image?', ...
              'Next', 'Back', 'Skip this', 'Goto...', 'Random', 'Quit');
switch choice
    case 1  % Next
        i = min(numel(L), i+1);
    case 2  % Back
        i = max(1, i-1);
    case 3  % Skip this
        i = min(numel(L), i+1);
    case 4  % Goto...
        answ = inputdlg(sprintf('Go to index (1..%d):', numel(L)), 'Goto', 1, {num2str(i)});
        if ~isempty(answ)
            v = str2double(answ{1});
            if ~isnan(v)
                i = min(max(1, round(v)), numel(L));
            end
        end
    case 5  % Random
        i = randi(numel(L));
    case 6  % Quit
        break;
end


    


end

writeRows(outCsv, rows);
fprintf('Saved labels to %s\n', outCsv);
end

function writeRows(outCsv, rows)
if isempty(rows)
    T = cell2table(cell(0,5), 'VariableNames',{'image','x','y','w','h'});
else
    T = cell2table(rows, 'VariableNames',{'image','x','y','w','h'});
end
if ~exist(fileparts(outCsv),'dir'), mkdir(fileparts(outCsv)); end
writetable(T, outCsv);
end

function b = basename(p)
[~,b,~] = fileparts(p);
end

function viewer(outDir)
% VIEWER  Arrow-key viewer for Original / Overlay / Mask with quick links.
% Usage: viewer('results')  or  viewer('results_flake')

if nargin==0, outDir = ""; end
if outDir=="" || ~isfolder(outDir)
    % auto-pick newest results* folder
    d = dir(fullfile(pwd,'results*')); d = d([d.isdir]);
    assert(~isempty(d),'No results* folders found in %s', pwd);
    [~,ix] = max([d.datenum]); outDir = fullfile(d(ix).folder, d(ix).name);
    fprintf('Auto-using results folder: %s\n', outDir);
end
assert(isfolder(outDir), 'Results folder not found: %s', outDir);

% collect overlays
Dov = dir(fullfile(outDir,'overlays','*_overlay.png'));
assert(~isempty(Dov), 'No overlays in %s/overlays', outDir);

% table of paths
N = numel(Dov);
T = table('Size',[N 7], ...
    'VariableTypes', ["string","string","string","string","double","double","double"], ...
    'VariableNames', {'image','overlay','mask','json','n_flakes','mean_conf','median_area_um2'});

for k = 1:N
    base = erase(string(Dov(k).name), "_overlay.png");

    T.overlay(k) = string(fullfile(Dov(k).folder, Dov(k).name));
    T.mask(k)    = string(fullfile(outDir,'masks', base + "_mask.png"));
    T.json(k)    = string(fullfile(outDir,'json',  base + "_flakes.json"));
    T.image(k)   = "";             % fill from JSON if available
    T.n_flakes(k) = NaN; T.mean_conf(k) = NaN; T.median_area_um2(k) = NaN;

    if isfile(T.json(k))
        S = jsondecode(fileread(T.json(k)));

        % try to get original image path
        img1 = getStrField(S, {'summary','image'});
        img2 = getStrField(S, {'summary','image_rel_path'});
        if strlength(img1)>0 && isfile(img1)
            T.image(k) = img1;
        elseif strlength(img2)>0
            cand = fullfile(pwd, img2);
            if isfile(cand), T.image(k) = string(cand); end
        end

        % robust numbers (accept scalar, vector, or numeric-string)
        nf = getNumField(S, {'summary','n_flakes'});
        mc = getNumField(S, {'summary','mean_conf'});
        ma = getNumField(S, {'summary','median_area_um2'});

        if ~isnan(nf), T.n_flakes(k)       = nf; end
        if ~isnan(mc), T.mean_conf(k)      = mc; end
        if ~isnan(ma), T.median_area_um2(k)= ma; end

        % if counts missing, compute from flakes
        if isnan(T.n_flakes(k)) && isfield(S,'flakes')
            T.n_flakes(k) = numel(S.flakes);
        end
        if isnan(T.mean_conf(k)) && isfield(S,'flakes') && ~isempty(S.flakes)
            try T.mean_conf(k) = mean([S.flakes.confidence]); catch, end
        end
        if isnan(T.median_area_um2(k)) && isfield(S,'flakes') && ~isempty(S.flakes)
            try T.median_area_um2(k) = median([S.flakes.area_um2]); catch, end
        end
    end
end

% viewer state
state.i = 1; state.mode = 2;  % 1=Orig 2=Overlay 3=Mask
h = figure('Name','flake viewer','KeyPressFcn',@onKey, 'NumberTitle','off'); %#ok<NASGU>
draw();

    function draw()
        clf;
        switch state.mode
            case 1
                if strlength(T.image(state.i))>0 && isfile(T.image(state.i))
                    im = imread(T.image(state.i)); ttl = 'Original';
                else
                    im = imread(T.overlay(state.i)); ttl = 'Original (fallback=overlay)';
                end
            case 2
                im = imread(T.overlay(state.i)); ttl = 'Overlay';
            case 3
                if isfile(T.mask(state.i))
                    mk = imread(T.mask(state.i)); if size(mk,3)>1, mk = rgb2gray(mk); end
                    im = uint8(mk); ttl = 'Mask';
                else
                    im = imread(T.overlay(state.i)); ttl = 'Mask (fallback=overlay)';
                end
        end
        imshow(im);
        [~,b] = fileparts(char(T.overlay(state.i)));
        b = erase(b, "_overlay");
        title(sprintf('[%d/%d] %s  |  %s  |  flakes=%s  conf=%s  medA=%s µm^2', ...
            state.i, height(T), b, ttl, showNum(T.n_flakes(state.i),0), ...
            showNum(T.mean_conf(state.i),2), showNum(T.median_area_um2(state.i),1)));
        drawnow;
    end

    function onKey(~,ev)
        switch ev.Key
            case {'rightarrow','pagedown','space'}
                state.i = min(height(T), state.i+1); draw();
            case {'leftarrow','pageup','backspace'}
                state.i = max(1, state.i-1); draw();
            case 'o', state.mode = 1; draw();
            case 'v', state.mode = 2; draw();
            case 'm', state.mode = 3; draw();
            case 'j', if isfile(T.json(state.i)), edit(T.json(state.i)); end
            case 'c', winopen(fullfile(outDir,'crops'));
            case 'r', close(gcf); viewer(outDir);
            case 'g'   % Go to index
    prompt = sprintf('Go to image [1..%d]: ', height(T));
    v = str2double(inputdlg(prompt, 'Go to', 1, {num2str(state.i)}));
    if ~isnan(v)
        state.i = min(max(1, round(v)), height(T)); draw();
    end
case 'h'   % show help
    helpmsg = sprintf([ ...
        'Controls:\n' ...
        '  → / Space : next image\n' ...
        '  ← / Back  : previous image\n' ...
        '  Home/End  : first/last\n' ...
        '  o / v / m : Original / overlay / mask\n' ...
        '  j         : open JSON\n' ...
        '  c         : open crops folder\n' ...
        '  g         : go to index\n' ...
        '  r         : reload viewer\n' ...
        'Click the image window first to give it keyboard focus.']);
    helpdlg(helpmsg, 'Viewer help');

        end
    end
end

% ---------- helpers ----------
function s = getStrField(S, path)
s = "";
try
    for i=1:numel(path)
        if isfield(S, path{i}), S = S.(path{i}); else, s = ""; return; end
    end
    if isstring(S) || ischar(S), s = string(S); end
catch, s = ""; end
end

function x = getNumField(S, path)
x = NaN;
try
    for i=1:numel(path)
        if isfield(S, path{i}), S = S.(path{i}); else, return; end
    end
    if isnumeric(S)
        if isempty(S), x = NaN; else, x = double(S(1)); end  % first element if vector
    elseif isstring(S) || ischar(S)
        t = str2double(string(S));
        if ~isnan(t), x = t; end
    end
catch, x = NaN; end
end

function s = showNum(x, nd)
if isnan(x), s="n/a"; else, s = string(num2str(x, ['%0.' num2str(nd) 'f'])); end
end

function meta = readImageMeta(folder, base)
%READIMAGEMETA  Load sidecar meta JSON for an image (best-effort).
% Looks for base_meta.json or base.json (case-insensitive).
% Normalizes common fields to: image_id, chip_id, magnification, stage_x_mm, stage_y_mm.

cands = { sprintf('%s_meta.json', base), sprintf('%s.json', base) };
d = dir(folder);
names = {d.name};
meta = struct();

for i = 1:numel(cands)
    idx = find(strcmpi(names, cands{i}), 1);
    if ~isempty(idx)
        try
            raw = fileread(fullfile(folder, names{idx}));
            meta = jsondecode(raw);
            meta = normalizeMeta(meta);
            return;
        catch ME
            warning('readImageMeta:decodeFailed', ...
                'Failed to parse %s: %s', names{idx}, ME.message);
            meta = struct(); % continue trying other candidates
        end
    end
end
% not found → return empty struct
end

function m2 = normalizeMeta(m)
% Map a wide variety of possible keys → normalized schema.
m2 = struct();

% IDs
m2.image_id = pick(m, ["image_id","tile_id","name","basename"]);
m2.chip_id  = pick(m, ["chip_id","sample_id","wafer_id","substrate_id"]);

% Magnification can come as number or string "20x"
mag = pick(m, ["magnification","mag","objective","objective_mag"]);
if ischar(mag) || isstring(mag)
    tok = regexp(string(mag), '(\d+(\.\d+)?)\s*x', 'tokens', 'once');
    if ~isempty(tok), mag = str2double(tok{1}{1}); else, mag = []; end
end
if ~isempty(mag) && isnumeric(mag), m2.magnification = mag; end

% Stage positions (prefer mm)
sx = pick(m, ["stage_x_mm","stage_x","x_mm","stage_x_um","x_um"]);
sy = pick(m, ["stage_y_mm","stage_y","y_mm","stage_y_um","y_um"]);

% unit inference
if ~isempty(sx) && ~isempty(sy)
    if contains(lower(join(string(fieldnames(m))')),'_um')
        % likely in micrometers
        sx = double(sx) / 1000;
        sy = double(sy) / 1000;
    else
        sx = double(sx);
        sy = double(sy);
    end
    m2.stage_x_mm = sx;
    m2.stage_y_mm = sy;
end

% fallback aliases
if ~isfield(m2,'image_id') && isfield(m2,'chip_id')
    m2.image_id = m2.chip_id + "_" + string(now*1e5);
end
end

function v = pick(s, keys)
v = [];
for i = 1:numel(keys)
    k = char(keys(i));
    if isfield(s,k) && ~isempty(s.(k))
        v = s.(k);
        return;
    end
end
end

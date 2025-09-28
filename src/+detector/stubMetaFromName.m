function m = stubMetaFromName(base)
% Infer minimal meta from filename when no sidecar JSON exists.
m = struct();

% chip id (e.g., chipA or chip-001)
tok = regexp(base, '(chip[-_]?\w+)', 'tokens', 'once');
if ~isempty(tok), m.chip_id = tok{1}; end

% tile id (e.g., tile_001 or tile-12)
tok = regexp(base, '(tile[-_]?\d+)', 'tokens', 'once');
if ~isempty(tok), m.image_id = tok{1}; end

% magnification: "_x20", "-x50", or "20x"
tok = regexp(base, '[_-]x?(\d+)\b|(\d+)x\b', 'tokens', 'once');
if ~isempty(tok)
    cand = tok(~cellfun('isempty', tok));
    m.magnification = str2double(cand{1});
end
end

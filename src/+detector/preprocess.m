function out = preprocess(I, varargin)
%PREPROCESS Normalize image and prepare features (toolbox-friendly & toolbox-free).
% Accepts EITHER:
%   preprocess(I, struct('UsedChannels','RG', ...))
% OR
%   preprocess(I, 'UsedChannels','RG', ...)

% ----- defaults -----
opts = struct( ...
    'UsedChannels','RGB', ...
    'BackgroundCeiling',[], ...
    'PerChannelPercentiles',[1 99], ...
    'Flatfield',[] );

% ----- parse varargin into opts (supports struct OR name-value) -----
if ~isempty(varargin)
    if isstruct(varargin{1})
        user = varargin{1};
        f = fieldnames(user);
        for k=1:numel(f), opts.(f{k}) = user.(f{k}); end
    else
        % name-value pairs
        if mod(numel(varargin),2)~=0
            error('Name-value inputs must come in pairs.');
        end
        for k=1:2:numel(varargin)
            name = varargin{k}; val = varargin{k+1};
            if ~ischar(name) && ~isstring(name), error('Option names must be char/string.'); end
            name = char(name);
            if ~isfield(opts, name)
                % be permissive: add unknown names to opts
                opts.(name) = val;
            else
                opts.(name) = val;
            end
        end
    end
end

% ----- ensure 3-channel uint8 -----
mustBeNumeric(I);
if isempty(I), error('Input image I is empty.'); end
if ndims(I)==2, I = repmat(I,1,1,3); end
I = I(:,:,1:3);
if ~isa(I,'uint8'), I = uint8(max(0,min(255,round(double(I))))); end
[H,W,~] = size(I);

% ----- optional ceiling -----
if ~isempty(opts.BackgroundCeiling)
    I = min(I, uint8(max(0, min(255, round(opts.BackgroundCeiling)))));
end

% ----- optional flatfield (multiplicative) -----
if ~isempty(opts.Flatfield)
    F = single(opts.Flatfield);
    if size(F,1)~=H || size(F,2)~=W, error('Flatfield size mismatch.'); end
    if size(F,3)==1, F = repmat(F,1,1,3); end
    I = uint8( max(0,min(255, round( 255 * ( single(I)/255 ./ max(F, eps('single')) ) ))) );
end

% ----- per-channel percentile stretch -----
pc = max(0,min(100,opts.PerChannelPercentiles));
Iproc = I;
for c=1:3
    v = single(I(:,:,c));
    if exist('prctile','file')==2
        lo = prctile(v(:), pc(1));
        hi = prctile(v(:), pc(2));
    else
        lo = percentile_(v(:), pc(1));
        hi = percentile_(v(:), pc(2));
    end
    if isfinite(lo) && isfinite(hi) && hi>lo
        vv = (v - lo) ./ (hi - lo);
        vv = max(0,min(1,vv));
        Iproc(:,:,c) = uint8(round(255*vv));
    end
end

% ----- channel subset -----
map = struct('R',1,'G',2,'B',3);
chs = upper(char(opts.UsedChannels));
chs = chs(ismember(chs,'RGB'));
if isempty(chs), chs = 'RGB'; end
usedIdx = arrayfun(@(ch) map.(char(ch)), num2cell(chs));
Iproc = Iproc(:,:,usedIdx);

% ----- outputs -----
out.Iproc   = Iproc;                                    % uint8 HxWxC
out.X       = single(reshape(Iproc,[],numel(usedIdx))) / 255; % N x C in [0,1]
out.usedIdx = usedIdx(:)';
out.usedStr = char(chs);
end

% ------- local fallback percentile (toolbox-free) -------
function q = percentile_(x, p)
x = double(x(:));
n = numel(x);
if n==0 || ~isfinite(p), q = NaN; return; end
p = max(0,min(100,double(p)));
xs = sort(x);
if n==1, q = xs; return; end
pos = 1 + (n-1)*(p/100); lo = floor(pos); hi = ceil(pos); w = pos - lo;
lo = max(1,min(n,lo)); hi = max(1,min(n,hi));
q = (1-w)*xs(lo) + w*xs(hi);
end

function R = evaluate_detections(resultsDir, gtCsv, iouThresh)
% EVALUATE_DETECTIONS  Compare detector JSON vs GT boxes (IoU).
% Only evaluates images present in GT (unlabeled images are skipped).
% Returns a clean numeric per-image table.
% Usage:
%   R = tools.evaluate_detections('results', fullfile('gt','gt_boxes.csv'), 0.30)

if nargin<3, iouThresh = 0.30; end
assert(isfolder(resultsDir), 'Results folder not found: %s', resultsDir);
assert(isfile(gtCsv), 'GT CSV not found: %s', gtCsv);

GT = readtable(gtCsv);

% ---- normalize GT image names to basenames (e.g., "10_x20")
GTbase = arrayfun(@(s) normBase(string(s)), GT.image);
GTset  = unique(GTbase);

Jdir = fullfile(resultsDir,'json');
D = dir(fullfile(Jdir,'*_flakes.json'));
assert(~isempty(D),'No JSON detections found in %s', Jdir);

tp=0; fp=0; fn=0;
per = struct('image',{}, 'tp',{}, 'fp',{}, 'fn',{});  % empty struct array

for k = 1:numel(D)
    % --- read detections
    S = jsondecode(fileread(fullfile(D(k).folder,D(k).name)));

    % --- choose best image identifier from JSON, then normalize
    imgKey = "";
    if isfield(S,'summary')
        if isfield(S.summary,'image_rel_path') && ~isempty(S.summary.image_rel_path)
            imgKey = string(S.summary.image_rel_path);
        elseif isfield(S.summary,'image') && ~isempty(S.summary.image)
            imgKey = string(S.summary.image);
        elseif isfield(S.summary,'image_base') && ~isempty(S.summary.image_base)
            imgKey = string(S.summary.image_base);
        end
    end
    if imgKey == "", imgKey = string(D(k).name); end
    jb = normBase(imgKey);  % e.g., "10_x20"

    % --- skip if no GT for this image
    if ~any(GTset == jb), continue; end

    % --- GT rows for this image
    gt = GT(GTbase == jb, :);

    % --- detections -> boxes [x y w h]
    detB = zeros(0,4);
    if isfield(S,'flakes') && ~isempty(S.flakes)
        db = vertcat(S.flakes.bbox);
        if ~isempty(db)
            detB = [[db.x].' [db.y].' [db.w].' [db.h].'];
        end
    end

    used = false(size(detB,1),1); localTP=0; localFN=0;
    for i=1:height(gt)
        g = [gt.x(i) gt.y(i) gt.w(i) gt.h(i)];
        ious = bboxIoU(detB, g);
        if isempty(ious), localFN = localFN + 1; continue; end
        [m, j] = max(ious);
        if m >= iouThresh && ~used(j)
            localTP = localTP + 1; used(j) = true;
        else
            localFN = localFN + 1;
        end
    end
    localFP = sum(~used);

    tp = tp + localTP; fp = fp + localFP; fn = fn + localFN;
    per(end+1) = struct('image', jb + ".png", 'tp', localTP, 'fp', localFP, 'fn', localFN); %#ok<AGROW>
end

prec = tp / max(tp+fp,1);
rec  = tp / max(tp+fn,1);
F1   = 2*prec*rec / max(prec+rec,eps);
fprintf('IoU>=%.2f  (evaluated %d labeled images)  TP=%d  FP=%d  FN=%d  →  P=%.3f  R=%.3f  F1=%.3f\n', ...
    iouThresh, numel(per), tp, fp, fn, prec, rec, F1);

if isempty(per)
    perTbl = cell2table(cell(0,4), 'VariableNames', {'image','tp','fp','fn'});
else
    perTbl = struct2table(per);
end
R = struct('tp',tp,'fp',fp,'fn',fn,'precision',prec,'recall',rec,'F1',F1,'iou',iouThresh,'perImage',perTbl);
end

% ---------- helpers ----------
function s = normBase(pathOrName)
% return lowercase base like "10_x20" with common suffixes removed
p = char(string(pathOrName));
[~,nm,~] = fileparts(p);
nm = lower(nm);
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

function iou = bboxIoU(B, g)
if isempty(B), iou = zeros(0,1); return; end
ix1 = max(B(:,1), g(1)); iy1 = max(B(:,2), g(2));
ix2 = min(B(:,1)+B(:,3), g(1)+g(3));
iy2 = min(B(:,2)+B(:,4), g(2)+g(4));
iw = max(0, ix2-ix1); ih = max(0, iy2-iy1);
inter = iw.*ih;
areaB = B(:,3).*B(:,4);
areag = g(3)*g(4);
iou = inter ./ max(areaB + areag - inter, eps);
end

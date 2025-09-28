function files = natsortfiles(files)
%NATSORTFILES natural sort for file paths (numbers in names sorted numerically).
% Minimal dependency-free implementation.

[pth,nam,ext] = cellfun(@fileparts, files, 'uni',0);
tokens = cellfun(@(s) regexp(s,'\d+|\D+','match'), nam, 'uni',0);

% pad tokens to same length
maxT = max(cellfun(@numel,tokens));
tokPad = cellfun(@(t) [t, repmat({''},1,maxT-numel(t))], tokens, 'uni',0);
% build sort keys: numeric tokens -> double, text -> string (lower)
Knum = zeros(numel(files), maxT);
Ktxt = strings(numel(files), maxT);
isNumTok = false(maxT,1);

for j=1:maxT
    col = cellfun(@(t) t{j}, tokPad, 'uni',0);
    isNum = cellfun(@(x) ~isempty(regexp(x,'^\d+$','once')), col);
    isNumTok(j) = any(isNum);
    Knum(:,j) = 0;
    Ktxt(:,j) = "";
    if any(isNum)
        Knum(isNum,j) = cellfun(@str2double, col(isNum));
        Ktxt(~isNum,j) = string(lower([col{~isNum}]));
    else
        Ktxt(:,j) = string(lower([col{:}]));
    end
end

% build a table of sort columns alternating numeric/text
S = [];
for j=1:maxT
    if isNumTok(j)
        S = [S, num2cell(Knum(:,j))]; %#ok<AGROW>
    else
        S = [S, cellstr(Ktxt(:,j))]; %#ok<AGROW>
    end
end

[~,ord] = sortrows(S);
files = files(ord);
end

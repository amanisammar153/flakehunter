function txt = jsonencodePretty(S)
%JSONENCODEPRETTY Pretty JSON using jsonencode + simple formatting.
raw = jsonencode(S);
% Insert newlines after commas/braces (simple, robust for our payload sizes)
txt = regexprep(raw, '(\{)', "$1\n");
txt = regexprep(txt, '(\})', "\n$1");
txt = regexprep(txt, ',"', ",\n""");
% Indent based on braces
lines = splitlines(txt);
indent = 0; sp = '  ';
for i=1:numel(lines)
    L = strtrim(lines{i});
    if startsWith(L, '}'), indent = max(0, indent-1); end
    lines{i} = [repmat(sp,1,indent) L];
    if endsWith(L, '{'), indent = indent + 1; end
end
txt = strjoin(lines, newline);
end

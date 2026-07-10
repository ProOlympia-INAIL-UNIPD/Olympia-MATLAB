function renameKeysInJSON(inputFile, oldStr, newStr, outputFile)
% RENAMEKEYSINJSON Replace substring inside JSON keys only (robust)
%
%   renameKeysInJSON(inputFile, oldStr, newStr)
%   renameKeysInJSON(inputFile, oldStr, newStr, outputFile)
%
% Replaces oldStr with newStr only inside JSON object keys (the "..." that
% are immediately followed by optional whitespace and a colon).
% The routine parses the JSON text char-by-char and handles escaped quotes.
%
% Limitations: molto robusta per JSON "standard". Non tenta di validare il
% JSON ma riconosce correttamente le stringhe e i backslash-escape.

    if nargin < 4
        outputFile = inputFile;
    end
    if ~isfile(inputFile)
        error('File "%s" not found.', inputFile);
    end

    txt = fileread(inputFile);
    n = length(txt);
    i = 1;
    out = strings(0);  % accumulate in string array pieces for speed
    pieceStart = 1;

    while i <= n
        ch = txt(i);
        if ch == '"'
            % append chunk before this string
            if pieceStart <= i-1
                out(end+1) = txt(pieceStart:i-1); %#ok<AGROW>
            end

            % parse the JSON string starting at i
            j = i + 1;
            unescaped = ""; % unescaped content
            rawSegmentStart = i; % original start for fallback
            while j <= n
                c = txt(j);
                if c == '\'
                    % escape sequence: take next char (if exists)
                    if j == n
                        % malformed escape at EOF, treat backslash literally
                        unescaped = unescaped + '\';
                        j = j + 1;
                    else
                        esc = txt(j+1);
                        % convert common escapes to their actual char
                        switch esc
                            case '"'
                                unescaped = unescaped + '"';
                            case '\'
                                unescaped = unescaped + '\';
                            case '/'
                                unescaped = unescaped + '/';
                            case 'b'
                                unescaped = unescaped + char(8);
                            case 'f'
                                unescaped = unescaped + char(12);
                            case 'n'
                                unescaped = unescaped + char(10); %newline
                            case 'r'
                                unescaped = unescaped + char(13);
                            case 't'
                                unescaped = unescaped + char(9);
                            case 'u'
                                % unicode escape \uXXXX : keep as-is sequence (approx)
                                if j+5 <= n
                                    hex = txt(j+2:j+5);
                                    % try to convert if hex valid
                                    try
                                        code = sscanf(hex, '%4x');
                                        if ~isempty(code)
                                            unescaped = unescaped + char(code);
                                        else
                                            % fallback: keep raw sequence
                                            unescaped = unescaped + ['\u' hex];
                                        end
                                    catch
                                        unescaped = unescaped + ['\u' hex];
                                    end
                                    j = j + 4; % additional increment below will add 1 more
                                end
                            otherwise
                                % unknown escape, keep the escaped char
                                unescaped = unescaped + esc;
                        end
                        j = j + 2;
                        continue;
                    end
                elseif c == '"'
                    % end of string
                    break;
                else
                    unescaped = unescaped + c;
                    j = j + 1;
                end
            end

            if j > n
                % unterminated string: fallback to copy rest and break
                out(end+1) = txt(i:end); %#ok<AGROW>
                pieceStart = n+1;
                break;
            end

            % j points at closing quote
            closingQuotePos = j;
            % lookahead after closing quote to see if this is a key (skip spaces) 
            k = closingQuotePos + 1;
            while k <= n && isspace(txt(k))
                k = k + 1;
            end
            isKey = (k <= n && txt(k) == ':');

            if isKey
                % do the replacement on the unescaped content
                newUnescaped = strrep(char(unescaped), oldStr, newStr);
                % re-escape backslashes and quotes and control chars to produce JSON-safe string
                escaped = matlabEscapeJSONString(newUnescaped);
                out(end+1) = ['"' escaped '"']; %#ok<AGROW>
            else
                % not a key: append the original substring from i to closingQuotePos
                out(end+1) = txt(i:closingQuotePos); %#ok<AGROW>
            end

            % advance pointers
            i = closingQuotePos + 1;
            pieceStart = i;
        else
            i = i + 1;
        end
    end

    % append any remaining piece
    if pieceStart <= n
        out(end+1) = txt(pieceStart:end);
    end

    newTxt = strjoin(out, '');

    % write to file
    fid = fopen(outputFile, 'w');
    if fid == -1
        error('Cannot open "%s" for writing.', outputFile);
    end
    fwrite(fid, newTxt, 'char');
    fclose(fid);

    fprintf('Keys replaced and saved to "%s"\n', outputFile);

    function s = matlabEscapeJSONString(s)
        % Escapes a matlab string for JSON string literal
        % - replaces backslash and " and control chars with JSON escapes
        s = char(s); % ensure char array
        % replace backslash first
        s = strrep(s, '\', '\\');
        s = strrep(s, '"', '\"');
        % replace control characters
        s = strrep(s, sprintf('\b'), '\b');
        s = strrep(s, sprintf('\f'), '\f');
        s = strrep(s, sprintf('\n'), '\n');
        s = strrep(s, sprintf('\r'), '\r');
        s = strrep(s, sprintf('\t'), '\t');
        % other non-printable chars could be encoded with \uXXXX if needed
    end

end


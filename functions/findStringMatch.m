function index = findStringMatch(cellArray, str1, str2)
% La funzione restituisce l'indice della cella che contiene entrambe le stringhe str1 e str2
index = [];
if iscell(str2)
    str2 = str2{:};
end
for i = 1:length(cellArray)
    % Verifica se entrambe le stringhe sono contenute nell'elemento della cella
    if contains(cellArray{i}, str1) && contains(cellArray{i}, str2)
        index = i; % Restituisci l'indice
    end
end
end
function newStruct = removeFieldSuffix(inputStruct, suffix)
    % removeFieldSuffix: Rimuove un suffisso dai nomi dei campi di una struttura, inclusi livelli nidificati e array di strutture.
    % INPUT:
    %   inputStruct - Struttura di ingresso con campi di nome variabile
    %   suffix      - Suffisso da rimuovere (es. '_att')
    % OUTPUT:
    %   newStruct   - Copia della struttura con i nomi dei campi modificati

    % Verifica che il suffisso sia una stringa
    if ~ischar(suffix) && ~isstring(suffix)
        error('Suffix must be a string.');
    end

    % Inizializza la nuova struttura vuota
    if numel(inputStruct) > 1 && isstruct(inputStruct)
        % Gestione degli array di strutture
        newStruct = repmat(struct(), size(inputStruct));
        for idx = 1:numel(inputStruct)
            newStruct(idx) = removeFieldSuffix(inputStruct(idx), suffix);
        end
    else
        % Ottieni i nomi dei campi della struttura di ingresso
        fieldNames = fieldnames(inputStruct);

        % Itera su ciascun campo
        newStruct = struct(); % Nuova struttura
        for i = 1:length(fieldNames)
            % Nome del campo corrente
            oldField = fieldNames{i};

            % Rimuovi il suffisso specificato se presente
            newField = regexprep(oldField, [suffix '$'], '');

            % Ottieni il valore del campo corrente
            fieldValue = inputStruct.(oldField);

            % Verifica se il valore è una struttura (per applicare la ricorsione)
            if isstruct(fieldValue)
                % Ricorsione per i campi nidificati
                newStruct.(newField) = removeFieldSuffix(fieldValue, suffix);
            else
                % Copia il valore normalmente
                newStruct.(newField) = fieldValue;
            end
        end
    end
end
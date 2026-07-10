function plot_with_text(x,y, value, intervallo,colore)

    hold on;
    % Calcolo del numero totale di punti
    num_punti = length(y);

    % Ciclo per posizionare il testo ogni intervallo punti
    for i = 1:intervallo:num_punti
        % Verifica se il vettore text contiene elementi
        if ~isempty(value)
            % Posizionamento del testo
            text(x(i), y(i), value, 'FontSize', 4,'BackgroundColor', 'white','HorizontalAlignment', 'center','Color',colore,'Margin',2);
        end
    end
end
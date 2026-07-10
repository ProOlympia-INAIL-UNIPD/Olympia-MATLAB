function letters = alphabet(N)
            % Check for N <= A-Z interval (65-90)
            if N > 26
                error('N must be <= 26');
            end

            % Generate letters from 'A' to 'Z'
            letters = arrayfun(@(x) string(char(x)), 65:(64+N));
        end
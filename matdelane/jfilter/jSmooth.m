function y = jSmooth(x, w)

    y = zeros(size(x));
    d = ndims(x);

    if d == 2

        bl = mean(x(:, 1:500), 2);
        n1 = size(x, 1);

        parfor ik = 1:n1

            y(ik, :) = bl(ik) + detrend(smooth(x(ik, :), w));

        end

    elseif d == 3
        
        bl = mean(x(:, :, 1:500), 3);
        n1 = size(x, 1);
        n2 = size(x, 2);

        parfor ik = 1:n1

            for jk = 1:n2
            
                y(ik, jk, :) = bl(ik, jk) + detrend(smooth(x(ik, jk, :), w));

            end

        end

    end

end
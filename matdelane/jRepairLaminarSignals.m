function y = jRepairLaminarSignals(x, method)

%{
y = jRepairLaminarSignals(x, cw, tw, method)

>Input args
x : Signal (<N x M> numeric matrix, N-channels with length M)
cw : Mean kernel rows (channel-widnow)
tw : Mean kernel columns (time-window)
method : Core operation function ("interp1")
 
>Output args
y : Repaired signal (<N x M> numeric matrix)

%}

    if ~exist('method', 'var')

        method = "interpolation";

    end

    if strcmpi(method, "interp")

        y = jInterpolation2D(x, cw, tw);

    elseif strcmpi(method, "median")

    end

    n = size(x, 1);
    m = size(x, 2);
    y = zeros(n, m);
    
    y(1, :) = x(1, :);
    y(n, :) = x(n, :);
    
    parfor i = 1:n

        if i <= cw

            for j = tw:m-tw
        
                y(i, j) = mean(x(i:i+cw, j-tw+1:j+tw-1), "all");
        
            end

        elseif i <= n-cw

            for j = tw:m-tw
        
                y(i, j) = mean(x(i-cw:i+cw, j-tw+1:j+tw-1), "all");
        
            end

        else

            for j = tw:m-tw
        
                y(i, j) = mean(x(i-cw:i, j-tw+1:j+tw-1), "all");
        
            end
            
        end

    end

end
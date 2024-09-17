function y = jMedianFilt2(x, cw, tw)

%{
y = jMedianFilt2(x, cw, tw)

>Input args
x : Signal (<N x M> numeric matrix)
cw : Median kernel rows (channel-widnow)
tw : Median kernel columns (time-window)
 
>Output args
y : Median-filtered signal (<N x M> numeric matrix)

%}

    n = size(x, 1);
    m = size(x, 2);
    y = zeros(n, m);
    
    y(1, :) = x(1, :);
    y(n, :) = x(n, :);
    
    for i = 1:n

        if i <= cw

            for j = tw:m-tw
        
                y(i, j) = median(x(i:i+cw, j-tw+1:j+tw-1), "all");
        
            end

        elseif i <= n-cw

            for j = tw:m-tw
        
                y(i, j) = median(x(i-cw:i+cw, j-tw+1:j+tw-1), "all");
        
            end

        else

            for j = tw:m-tw
        
                y(i, j) = median(x(i-cw:i, j-tw+1:j+tw-1), "all");
        
            end
            
        end

    end

end

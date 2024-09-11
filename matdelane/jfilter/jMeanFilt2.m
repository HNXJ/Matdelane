function y = jMeanFilt2(x, cw, tw)

%{
y = jMeanFilt2(x, cw, tw)

>Input args
x : Signal (<N x M> numeric matrix)
cw : Mean kernel rows (channel-widnow)
tw : Mean kernel columns (time-window)
 
>Output args
y : Filtered signal (<N x M> numeric matrix)

%}

    n = size(x, 1);
    m = size(x, 2);
    y = zeros(n, m);
    
    y(1, :) = x(1, :);
    y(n, :) = x(n, :);
    
    for i = 1:n

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
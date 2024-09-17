function y = jApproximateEntropy(x)

    N = size(x, 1);
    y = zeros(1, N);
    
    for i = 1:N
    
        y(i) = approximateEntropy(x(i, :));
    
    end

end
function d = couplingDistance(A, B)
% Computes mean nearest distance from peaks A to peaks B. B is the reference, so one nearest distance to each B peak is taken into account regardless of the number of A peaks. 
    
    dists = zeros(size(B));

    for i = 1:1:length(B)
        d2 = [];
        for j = 1:1:length(A)
            d2(j) = abs(A(j) - B(i));
        end     
        dists(i) = min(d2);
    end

    d = mean(dists);
end

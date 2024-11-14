function D = Dmatrix_Fourier(N, xjs)
    D = zeros(N, N);
    for i = 1:N
        for j = 1:N
            if i == j
                D(i, j) = 0;
            else
                D(i, j) = diff_h(N, j, xjs(i), xjs);
            end
        end
        s = sum(D(i, :));
        D(i, i) = -s;
    end
end
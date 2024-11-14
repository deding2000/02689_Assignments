function dh = diff_h(N, j, x, xjs)
    dx = x - xjs(j);
    dh = cos(N*dx/2)/tan(dx/2)/2;
    dh = dh + sin(N*dx/2)*(-(1 + 1/tan(dx/2)^2)/2)/N;
end
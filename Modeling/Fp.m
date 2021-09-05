function y = Fp(v,r)
    % air blow propulsion force
    if r == 0.01
       y = 0.0000019 * v^3 - 0.0000846 * v^2 + 0.0013006 * v - 0.0020981;
    else
        y = 0;
    end
end
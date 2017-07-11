function [AA, B, status] = solve_using_solab(G0, G1, n_pre, n_jump, n_shock)
    [ff, p, status] = solab(G0, G1, n_pre);
    if (status == -1)
        AA = [];
        B = [];
        return;
    end
    
    if (n_pre == 2) 
        B = [1; p(2,1)/p(1,1); ff(1)];
        AA = [p zeros(2,1); ff*p 0];
    elseif (n_pre == 3)
        AA = [p zeros(3,1); ff * p, 0];
        B = (1./p(1,1)) * [p(:,1); ff * p(:,1)];
    else
        AA = [p zeros(n_pre,n_jump); ff*p zeros(n_jump,n_jump)];
        B = zeros(n_pre+n_jump, n_shock);
        B(1:n_shock,1:n_shock) = eye(n_shock);
        for idx = n_shock+1:n_pre
            for shock_idx = 1:n_shock
                B(idx,shock_idx) = p(idx,shock_idx) / p(shock_idx,shock_idx);
            end
        end
        for shock_idx = 1:n_shock
            for idx = 1:n_jump
                B(n_pre+idx,shock_idx) = ff(idx,:) * p(:,shock_idx) ./ p(shock_idx,shock_idx);
            end
        end
    end
end
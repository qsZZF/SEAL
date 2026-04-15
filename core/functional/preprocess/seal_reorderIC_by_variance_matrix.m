function ica = seal_reorderIC_by_variance_matrix(ica)
            v = var(double(ica.S), 0, 2);
            [~, order] = sort(v, 'descend');
            ica.S = ica.S(order, :);
            ica.W = ica.W(order, :);
            ica.A = ica.A(:, order);
        end
%     Copyright (C) 2026, University of Oslo  [Mathias Hudoba de Badyn]
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.

function plot_weighted_triangulated_annulus(E, F, W0, W1, W2)
    % Plot a weighted triangulated annulus
    %
    % Inputs:
    %   E  - ne x 2 edge list
    %   F  - nf x 3 face list
    %   W0 - nv x nv diagonal matrix of node weights
    %   W1 - ne x ne diagonal matrix of edge weights
    %   W2 - nf x nf diagonal matrix of face weights
    %
    % The layout is computed from a spectral (Laplacian) drawing.

    nv = size(W0,1);
    ne = size(E,1);
    nf = size(F,1);

    w0 = diag(W0);
    w1 = diag(W1);
    w2 = diag(W2);

    %% Build unweighted graph adjacency for layout
    A = zeros(nv,nv);
    for k = 1:ne
        i = E(k,1);
        j = E(k,2);
        A(i,j) = 1;
        A(j,i) = 1;
    end

    %% Spectral / Laplacian layout
    deg = sum(A,2);
    L = diag(deg) - A;

    % Use the 2nd and 3rd eigenvectors for 2D coordinates
    [U, D] = eig(L);
    [~, idx] = sort(diag(D), 'ascend');

    x = U(:, idx(2));
    y = U(:, idx(3));

    % Normalize coordinates
    x = x - mean(x);
    y = y - mean(y);
    scale = max(max(abs(x)), max(abs(y)));
    if scale > 0
        x = x / scale;
        y = y / scale;
    end

    % Optional mild radial spreading for prettier separation
    r = sqrt(x.^2 + y.^2);
    rmax = max(r);
    if rmax > 0
        x = x ./ (r + 0.15);
        y = y ./ (r + 0.15);
    end

    %% Visual scaling
    node_sizes = rescale(w0, 80, 400);     % marker areas
    edge_widths = rescale(w1, 0.5, 6);     % line widths
    face_vals   = rescale(w2, 0.15, 1.0);  % color intensity

    %% Plot
    figure;
    hold on;
    axis equal off;

    % Draw faces first
    cmap = parula(256);
    for f = 1:nf
        verts = F(f,:);
        Xf = x(verts);
        Yf = y(verts);

        cidx = max(1, min(256, round(1 + (size(cmap,1)-1)*(face_vals(f)-min(face_vals)) / max(eps, max(face_vals)-min(face_vals)))));
        patch(Xf, Yf, cmap(cidx,:), ...
            'FaceAlpha', 0.65, ...
            'EdgeColor', 'none');
    end

    % Draw edges
    for e = 1:ne
        i = E(e,1);
        j = E(e,2);
        plot([x(i), x(j)], [y(i), y(j)], 'k-', ...
            'LineWidth', edge_widths(e));
    end

    % Draw nodes
    scatter(x, y, node_sizes, 'w', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

    % % Node labels
    % for i = 1:nv
    %     text(x(i), y(i), sprintf('  %d', i), ...
    %         'FontSize', 11, 'FontWeight', 'bold', ...
    %         'HorizontalAlignment', 'left', ...
    %         'VerticalAlignment', 'middle');
    % end

    colormap(cmap);
    cb = colorbar;
    cb.Label.String = 'Face weight';
    caxis([min(w2), max(w2)]);
    cb.Label.FontSize = 18;    % label font size
    cb.FontSize = 16;          % numbers on the colorbar


    hold off;
end
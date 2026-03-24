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

function [B1, B2, B3] = construct_boundary_matrices(V, E, F, T)
    % Constructs the boundary matrices B1 (vertex-edge), B2 (edge-face), and B3 (face-tetrahedron)
    % ensuring that B1 * B2 = 0 and B2 * B3 = 0 for a correctly oriented simplicial complex.
    %
    % Inputs:
    %   - V: Vertices (not used directly in construction but included for completeness)
    %   - E: Edges (each row is [v1, v2])
    %   - F: Faces (each row is [v1, v2, v3])
    %   - T: Tetrahedra (each row is [v1, v2, v3, v4])
    %
    % Outputs:
    %   - B1: Vertex-Edge Incidence Matrix
    %   - B2: Edge-Face Incidence Matrix
    %   - B3: Face-Tetrahedron Incidence Matrix

    nv = size(V, 1); % Number of vertices
    ne = size(E, 1); % Number of edges
    nf = size(F, 1); % Number of faces
    nt = size(T, 1); % Number of tetrahedra

    % Initialize B1 (Vertex-Edge Incidence Matrix)
    B1 = zeros(nv, ne);
    for i = 1:ne
        v1 = E(i,1);
        v2 = E(i,2);
        B1(v1, i) = -1; % Start vertex
        B1(v2, i) = 1;  % End vertex
    end

    % Initialize B2 (Edge-Face Incidence Matrix)
    B2 = zeros(ne, nf);
    for j = 1:nf
        % Extract vertices of the face
        v1 = F(j,1);
        v2 = F(j,2);
        v3 = F(j,3);
        
        % Define edges of the face in counterclockwise order
        face_edges = [v1 v2; v2 v3; v3 v1];

        % Assign edges to faces with correct orientation
        for k = 1:3
            edge = sort(face_edges(k,:)); % Sort to match stored edges
            edge_idx = find(ismember(E, edge, 'rows'));
            if ~isempty(edge_idx)
                if isequal(E(edge_idx,:), face_edges(k,:)) 
                    B2(edge_idx, j) = 1; % Correct orientation
                else
                    B2(edge_idx, j) = -1; % Flip orientation
                end
            end
        end
    end

    % Initialize B3 (Face-Tetrahedron Incidence Matrix)
    B3 = zeros(nf, nt);
    face_signs = [1, -1, 1, -1]; % Alternating signs to enforce outward normals

    for k = 1:nt
        % Extract vertices of the tetrahedron
        v1 = T(k,1);
        v2 = T(k,2);
        v3 = T(k,3);
        v4 = T(k,4);
        
        % Define the four faces of the tetrahedron
        faces = [v1 v2 v3; v1 v2 v4; v1 v3 v4; v2 v3 v4];

        % Assign faces to tetrahedra with correct orientation
        for f = 1:4
            face = sort(faces(f,:)); % Sort to match stored faces
            face_idx = find(ismember(F, face, 'rows'));
            if ~isempty(face_idx)
                B3(face_idx, k) = face_signs(f); % Assign with correct sign
            end
        end
    end

    % Verify B1 * B2 == 0
    error_B1B2 = norm(B1 * B2);
    if error_B1B2 > 1e-10
        warning('B1 * B2 is not zero! Possible orientation issues.');
    end

    % Verify B2 * B3 == 0
    error_B2B3 = norm(B2 * B3);
    if error_B2B3 > 1e-10
        warning('B2 * B3 is not zero! Possible orientation issues.');
    end
end

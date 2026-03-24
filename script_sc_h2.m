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

clear all; close all; clc;

addpath(genpath('./'))

% this script calculates the H2 norm of a simplicial complex, and verifies
% the notions of effective resistance

% s^1 = B_1^T s_0 + s_H^1 + B_2 s^2 is a signal defined on the edges

rng(1337)

%% Define the Simplicial Complex (Tetrahedral Mesh)

% simple tetrahedron

V = [ 0;
      1;
      0;
      0 ]; % 4 Vertices

E = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % 6 Edges
F = [1 2 3; 1 2 4; 1 3 4; 2 3 4]; % 4 Faces (Triangles)
T = [1 2 3 4]; % 1 Tetrahedron

% two tetrahedra joined by a face
% Vertices
V = [0; 1; 0; 0; 1];  % 5 Vertices (coordinates are placeholders)

% Edges (unordered pairs, assumed to be sorted uniquely elsewhere if needed)
E = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4;   % edges of first tetrahedron
     1 5; 2 5; 3 5];                 % edges for the second tetrahedron

% Faces (Triangles)
F = [1 2 3;   % shared face
     1 2 4;
     1 3 4;
     2 3 4;   % faces of first tetrahedron
     1 3 5;
     1 2 5;
     2 3 5];  % additional faces for second tetrahedron

% Tetrahedra
T = [1 2 3 4;
     1 2 3 5];  % two tetrahedra joined by face [1 2 3]


%%%%%%%%%%%% Triangulated annulus: b1 = 1

V = [0; 1; 2; 3; 4; 5; 6; 7];   % placeholder coordinates

% Outer square: 1-2-3-4
% Inner square: 5-6-7-8
E = [1 2;
     2 3;
     3 4;
     1 4;    % outer boundary

     5 6;
     6 7;
     7 8;
     5 8;    % inner boundary

     1 5;
     2 6;
     3 7;
     4 8;    % radial edges

     2 5;
     3 6;
     4 7;
     1 8];   % diagonals for triangulation

F = [1 2 5;
     2 5 6;
     2 3 6;
     3 6 7;
     3 4 7;
     4 7 8;
     1 4 8;
     1 5 8];

T = zeros(0,4);   % no tetrahedra

%%%%%%%%%%%% weights

nv = size(V,1);  % Number of vertices
ne = size(E,1);  % Number of edges
nf = size(F,1);  % Number of faces
nt = size(T,1);  % Number of tetrahedra

W0 = diag(0.5 + rand(nv,1));
W1 = diag(0.5 + rand(ne,1));
W2 = diag(0.5 + rand(nf,1));
W3 = diag(0.5 + rand(nt,1));

% W0 = eye(nv);
% W1 = eye(ne);
% W2 = eye(nf);
% W3 = eye(nt);

%% Step 2: Construct Boundary Matrices

[B1, B2, B3] = construct_boundary_matrices(V, E, F, T);

%check boundary maps
assert(norm(B1*B2)<1e-10)
assert(norm(B2*B3)<1e-10)

%% Compute Hodge Laplacians
L0 =                     inv(W0)*B1*W1*B1'; % Vertex Laplacian
L1 = B1'*inv(W0)*B1*W1 + inv(W1)*B2*W2*B2'; % Edge Laplacian
L2 = B2'*inv(W1)*B2*W2; % + inv(W2)*B3*W3*B3'; % Face Laplacian

L0up = inv(W0)*B1*W1*B1';
L1up = inv(W1)*B2*W2*B2';
L1dn = B1'*inv(W0)*B1*W1;
L2dn = B2'*inv(W1)*B2*W2;

%% Compute centering matrix

M1dn = B1'*pinv(L0up)*inv(W0)*B1*W1;
M1up = inv(W1)*B2*W2*pinv(L2dn)*B2';

M1 = B1'*pinv(L0up)*inv(W0)*B1*W1 + inv(W1)*B2*W2*pinv(L2dn)*B2';

%% Harmonic basis
H = null(L1);   % basis for ker(L1)

%% Weighted harmonic projector
PH = H * ((H' * W1 * H) \ (H' * W1));

%% Weighted disagreement projector
M2 = eye(ne) - PH;

%check agreement with the above centering matrix
assert(norm(M1-M2)<1e-12)

%% Two equivalent systems

% 1) centered disagreement dynamics
% xd_dot = -(L1 + PH) xd + Pi1 u
% y      = xd

As1 = -(L1 + PH);
Bs1 = M1;
Cs1 = sqrt(W1)*eye(ne);
Ds1 = zeros(ne);

sys = ss(As1,Bs1,Cs1,Ds1);
h2_sys1 = norm(sys,2)^2;

%% Separate the disagreement dynamics into up and down dynamics

% check
assert(norm(M1dn + M1up - M1) < 1e-10)

% split dynamics
A_dn = -(L1dn + PH + M1up);
B_dn = M1dn;
C_dn = sqrt(W1)*eye(ne);
D_dn = zeros(ne);

A_up = -(L1up + PH + M1dn);
B_up = M1up;
C_up = sqrt(W1)*eye(ne);
D_up = zeros(ne);

sys_dn = ss(A_dn,B_dn,C_dn,D_dn);
sys_up = ss(A_up,B_up,C_up,D_up);

h2_dn = norm(sys_dn,2)^2;
h2_up = norm(sys_up,2)^2;

assert(abs(norm(sys_dn,2)^2 + norm(sys_up,2)^2 - norm(sys,2)^2)<1E-10)

% up theoretical version
L1upsym = sqrtm(W1)*L1up/sqrtm(W1);
h2_sys1_up_theory = 0.5*trace(sqrtm(W1)*pinv(L1upsym)*sqrtm(W1));

assert(abs(norm(sys_up,2)^2 - h2_sys1_up_theory)<1E-10)

% dn theoretical version
L1dnsym = sqrtm(W1)*L1dn/sqrtm(W1);
h2_sys1_dn_theory = 0.5*trace(sqrtm(W1)*pinv(L1dnsym)*sqrtm(W1));

assert(abs(norm(sys_dn,2)^2 - h2_sys1_dn_theory)<1E-10)


%% effective resistance of a graph
% Useful square roots
S1  = sqrtm(W1);
S1i = inv(S1);

S2  = sqrtm(W2);
S2i = inv(S2);


H0 = null(L0);   
PH0 = H0*((H0'*W0*H0)\(H0'*W0));
As0 = -(L0 + PH0);
Bs0 = eye(nv) - PH0;
Cs0 = sqrtm(W0);

sys0 = ss(As0,Bs0,Cs0,0);
h20 = norm(sys0,2)^2;

%note that when node weights w0 = 0, then this is just the effective resistance.
L0up_sym = sqrtm(W0)*L0up/sqrtm(W0);
L0up_sym_pinv = pinv(L0up_sym);
Ginv = inv(L0up_sym + sqrtm(W0)*PH0/sqrtm(W0));
h20_theory = 0.5*trace(sqrtm(W0)*pinv(L0up_sym)*sqrtm(W0));
assert(abs(h20-h20_theory)<1E-12)

L0pinv = pinv(L0);
R0_true = zeros(nv,nv);

for ii = 1:nv
    for jj = 1:nv
        R0_true(ii,jj) = L0pinv(ii,ii) + L0pinv(jj,jj) - 2*L0pinv(ii,jj);
    end
end

% ignoring the weights from the nodes
assert(abs(0.5*ones(1,nv)*R0_true*ones(nv,1)/nv-trace(L0pinv))<1e-10 )

% with weights on nodes
trL0sym = nv*trace(L0up_sym_pinv);

R0_weight3 = zeros(nv,nv); 

for ii = 1:nv
    for jj = 1:nv
        
        eiminej = zeros(nv,1);
        if ii ~= jj
            eiminej(ii) = 1;
            eiminej(jj) = -1;
        end
        
        R0_weight3(ii,jj) = (sqrtm(W0)*eiminej)'*L0up_sym_pinv*(sqrtm(W0)*eiminej);
    end
end


assert(abs(0.5*sum(sum(R0_weight3))/(2*nv) - h20)<1E-12)

% try the Kirchoff index, doesn't work. Have to use the formula in R0_weight3
K=0;
for ii = 1:nv
    for jj = 1:ii-1
        eiminej = zeros(nv,1);
        eiminej(ii) = 1;
        eiminej(jj) = -1;
        
        K = K + eiminej'*W0*pinv(L0)*eiminej;
    end
end

%% effective resistance of faces on edge dynamics using up Laplacian

L1upsym_pinv = pinv(L1upsym);
h2_up_theory = 0.5*trace(sqrtm(W1)*L1upsym_pinv*sqrtm(W1));
assert(abs(h2_up-h2_up_theory)<1E-12)

R1_weight = zeros(ne,ne);

for ii = 1:ne
    for jj = 1:ne
        eiminej = zeros(ne,1);
        if ii ~= jj
            eiminej(ii) = 1;
            eiminej(jj) = -1;
        end
        R1_weight(ii,jj) = (sqrtm(W1)*eiminej)'*L1upsym_pinv*(sqrtm(W1)*eiminej);
    end
end

Gup = sqrtm(W1)*L1upsym_pinv*sqrtm(W1);

% what the Kirchoff index should equal
Ki_up1 = ne*trace(Gup) - ones(1,ne)*Gup*ones(ne,1);
assert(abs(h2_up-(Ki_up1 + ones(1,ne)*Gup*ones(ne,1))/(2*ne))<1E-12)

% Kirchoff index from R1_weight
Ki_up2 = 0.5*sum(sum(R1_weight));

assert(abs(Ki_up1-Ki_up2)<1E-12)

%% effective resistance of nodes on edge dynamics using dn Laplacian

L1dnsym_pinv = pinv(L1dnsym);
h2_dn_theory = 0.5*trace(sqrtm(W1)*L1dnsym_pinv*sqrtm(W1));
assert(abs(h2_dn-h2_dn_theory)<1E-12)

R1_dn_weight = zeros(ne,ne);

for ii = 1:ne
    for jj = 1:ne
        eiminej = zeros(ne,1);
        if ii ~= jj
            eiminej(ii) = 1;
            eiminej(jj) = -1;
        end
        R1_dn_weight(ii,jj) = (sqrtm(W1)*eiminej)'*L1dnsym_pinv*(sqrtm(W1)*eiminej);
    end
end

Gdn = sqrtm(W1)*L1dnsym_pinv*sqrtm(W1);

% what the Kirchoff index should equal
Ki_dn1 = ne*trace(Gdn) - ones(1,ne)*Gdn*ones(ne,1);
assert(abs(h2_dn-(Ki_dn1 + ones(1,ne)*Gdn*ones(ne,1))/(2*ne))<1E-12)

% Kirchoff index from R1_weight
Ki_dn2 = 0.5*sum(sum(R1_dn_weight));
assert(abs(Ki_dn1-Ki_dn2)<1E-12)
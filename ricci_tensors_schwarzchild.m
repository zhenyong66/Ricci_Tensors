% This program computes the Christoffel symbols of the second kind, Riemann
% tensors, and Ricci tensors for the Schwarzchild metric. It returns the
% non-zero Ricci tensors.

%% define matrix g of the metric and the variables
syms t r theta phi % create symbolic variables
syms a(r) b(r);
g = [-exp(2*a(r)) 0 0 0; 0 exp(2*b(r)) 0 0; 0 0 r^2 0; 0 0 0 r^2*(sin(theta))^2]; % matrix g for the Schwarzchild metric
G = inv(g);
X = [t r theta phi];

%% calculate and store christoffel symbols of 2nd kind
disp('Calculating Christoffel symbols of second kind...');
gamma_ij_l = 0;
for i = 1:4
    for j = 1:4
        for l = 1:4
            for k = 1:4
                Q = (1/2)*G(l,k)*(diff(g(j,k),X(i))+diff(g(i,k),X(j))-diff(g(i,j),X(k))); % Christoffel symbol of 1st kind
                gamma_ij_l = gamma_ij_l + Q;
            end
        CS2(i,j,l) = gamma_ij_l;
        gamma_ij_l = 0;
        end
    end
end

%% calculate and store Riemann curvature tensors
disp('Calculating Riemann curvature tensors...');
R_ijk_l = 0;
L_3 = 0;
for i = 1:4
    for j = 1:4
        for k = 1:4
            for l = 1:4
                L_1 = diff(CS2(i,k,l),X(j)) - diff(CS2(i,j,l),X(k));
                for m = 1:4
                    L_2 = CS2(i,k,m)*CS2(m,j,l)-CS2(i,j,m)*CS2(m,k,l);
                    L_3 = L_3 + L_2;
                end
                R_ijk_l = R_ijk_l + L_1 + L_3;
                Riemann(i,j,k,l) = R_ijk_l;
                L_3 = 0;
                R_ijk_l = 0;
            end
        end
    end
end

%% calculate and store Ricci tensors
disp('Calculating Ricci tensors...');
R_ij = 0;
for i = 1:4
    for j = 1:4
        for l=1:4
            R_ij = R_ij + Riemann(i,l,j,l);
        end
        Ricci(i,j) = R_ij;
        R_ij = 0;
    end
end

%% display results
disp('The non-zero Ricci tensors are:');
for i = 1:4
    for j = 1:4
        if Ricci(i,j) ~= 0;
            T = ['R_', num2str(i), num2str(j), ':'];
            disp(T);
            disp(Ricci(i,j));
        end
    end
end

disp('To call the following:');
disp('(i) Christoffel symbol of the second kind, type CS2(i,j,k),');
disp('(ii) Riemann tensor, type Riemann(i,j,k,l),');
disp('(iii) Ricci tensor, type Ricci(i,j).');
%%% This function is used to calculate the A_matrix and these  coefficiencies.
%%% the position of 1 in the right unit vector is depend on the
%%% value of m, besides, m can not be bigger than the sum of p and q,
%%% because of the limit the dimension of the matrix and the maximum 
%%% possible derivate these stencil points can approximate.

function [Coeff, A_matrix] = Coefficiency_fun(m_in, P_in, Q_in)
    size_matrix = P_in + Q_in + 1;          % the size of matrix is depend on the value of p and q. The exact relationship
                                            % between them are left.
    A_matrix = zeros(size_matrix);          % initialize the A_matrix
    Right_unit_vector = zeros(size_matrix, 1);          % initialize the unit vector
    Right_unit_vector(m_in + 1) = 1;                    % value the related element to 1
    
    for distance = -P_in :  Q_in            % solve the A_matrix with different distances and powers
        for power = 0 :  size_matrix - 1
            A_matrix(power + 1,distance + P_in + 1) = distance^power / factorial(power);
        end
    end
  
    Coeff = A_matrix \ Right_unit_vector;   % solve the Coefficiency vector
    return;
end
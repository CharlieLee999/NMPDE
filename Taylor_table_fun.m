%  With the input of m, p, q, to generate the Taylor table
function [Coeff_out, A_mat, Taylor_table_without_1, Taylor_table_with_1, Accuracy_r] = Taylor_table_fun(m_in, P_in, Q_in)
    size_matrix = P_in + Q_in + 1;              % The same with Coefficiency fun
    [Coeff_out, A_mat] = Coefficiency_fun(m_in,P_in,Q_in);              % Calling the Coefficiency fun to get Coeff and A_matrix  
    Taylor_table_without_1 = repmat(Coeff_out, 1,size_matrix).* A_mat'; % To generate the Taylor table without 1
    Right_unit_vector = zeros(1, size_matrix);                          % Initialize the unit vector
    Right_unit_vector(m_in + 1) = -1;                                   % The unit vector should be minus because it moved to the left of equa.
    Taylor_table_with_1 = [Right_unit_vector;Taylor_table_without_1];
        
    higher_column = zeros(1,size_matrix);                               % Row vector, To initialize and generate the higher Taylor expansion of p+q+1 order.                    
    for distance = -P_in :  Q_in
        higher_column(distance + P_in + 1) = distance^size_matrix / factorial(size_matrix);
    end
                                                                        % condition of equility 
    sum_higher_deri = sum(higher_column * Coeff_out);
    if abs(sum_higher_deri) <= 1e-5                                     % Whether the sum of higher taylor expansion is small enough, The criteria of convergence is 1e-5
        Accuracy_r = size_matrix +1 - m_in;                             % If yes, the accuracy should be higher
    else 
        Accuracy_r = size_matrix - m_in;                                % If no, the accuracy should be lower
    end 
   
    return;
end 


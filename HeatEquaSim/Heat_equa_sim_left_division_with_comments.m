%%% Receive input parameter: the size of mesh grid
N = input('Please input the positive integer N(size of mesh grid):');

t_0 = tic;
% Length = input('Please input the positive number L(length of space):')
% Time = input('Please input the positive number T(sim period):')

%%% Initialize mesh
a = 1;
Length = 1;
Time = 0.2;
tho = 1/N;
h = 1/N;
Space_step = Length / h;
Time_step = Time / tho;
Size_mat_mesh_space = Space_step + 1;
Size_mat_mesh_time = Time_step +1;


%%% Because tao == h and a == 1, lamda= N
Lamda  = a * (tho) / (h)^2;

%%% The size of meshing matrix is (Size_mat_mesh_space, Size_mat_mesh_time)
%%% Initialise and generate the Mesh_mat
Mesh_mat = zeros(Size_mat_mesh_time, Size_mat_mesh_space);
Space_vect = 0: h : Length;
Mesh_mat(1,:) = sin(pi * Space_vect);

%%% Generate and Calculate the coefficient matrix A and B 
Size_mat_A = Space_step - 1;
A_diag = (1 - Lamda) * ones(1,Size_mat_A);
A_diag_u_L = Lamda/2 * ones(1,Size_mat_A-1);
A_mat = diag(A_diag) + diag(A_diag_u_L,1) + diag(A_diag_u_L, -1);

% A_mat(1:Size_mat_A - 1, 2:Size_mat_A) = A_mat(1:Size_mat_A - 1, ...
%     2:Size_mat_A) + Lamda/2 * eye(Size_mat_A - 1);
% A_mat(2:Size_mat_A, 1:Size_mat_A - 1) = A_mat(2:Size_mat_A, ...
%     1:Size_mat_A - 1) + Lamda/2 * eye(Size_mat_A - 1);

B_mat = -A_mat + 2 * eye(Size_mat_A);



B_mat_inv = (B_mat)^-1;
%%% Calculate u(ih,nt)
for i = 2 : Size_mat_mesh_time
    Mesh_mat(i, 2: Space_step) = B_mat_inv * ( A_mat * Mesh_mat(i-1, 2:Space_step)');
end
fprintf('Running time:')
disp(toc(t_0))
%%% Make a video of the changing value of heat equation in different pooins
%%% and time.
video = VideoWriter('Heat_equ._sim.avi','Uncompressed AVI');
video.FrameRate = 10;       % Set by experience
open(video);

figure()
for i = 1: Size_mat_mesh_time
    title('Appr. and theo. heat equ. in diff. time')
    xlabel('Position')
    ylabel('Heat function')
    plot(Space_vect, Mesh_mat(i,:))
    hold on
    plot(Space_vect, sin(pi * Space_vect) * exp(- pi*pi*(i-1) * Time/ Time_step))
    hold on
    pause(0.05)
    legend(strcat('Time',num2str((i-1)*Time/Time_step)))
    F(i) = getframe(gcf);
end
saveas(gcf,'Appro_and_theo_values_in_diff_time','jpg')

writeVideo(video,F)
close(video)

%%% The max error when t equals 0.2
final_error = abs(Mesh_mat(Size_mat_mesh_time,:)-sin(pi * Space_vect) * exp(- pi*pi*Time));
final_max_error = max(final_error);
final_max_error_posi = find(final_error == final_max_error);

%%% To calculate the max error during the whole simulation period
error = zeros(1,Size_mat_mesh_time);
for i = 1: Size_mat_mesh_time
    error(i) = max(abs(Mesh_mat(i,:)-sin(pi * Space_vect) * exp(- pi*pi* ...
        (i-1) * Time/Time_step)));
end
max_error = max(error);
sim_step_num = find(error == max_error);
max_error_vector = abs(Mesh_mat(sim_step_num,:)-sin(pi * Space_vect)...
    * exp(- pi*pi* (sim_step_num-1) * Time/Time_step));
max_error_posi = find(max_error_vector == max_error);

%%% Print the result out
fprintf('The final max approximation error is:\t%f\n',final_max_error)
fprintf('It happens when x = %f \n\n',(final_max_error_posi-1)* Length/Space_step)

fprintf('The max approximation error during the whole simulation period is:\t %f\n', max_error)
fprintf('The simulation time t = %f \n', (sim_step_num-1)* Time/Time_step)
fprintf('It happens when x = %f \n\n',(max_error_posi-1)* Length/Space_step)

%%% Make a video of the changing value of approximation error in 
%%% different pooins and time
video2 = VideoWriter('App_error.avi','Uncompressed AVI');
video2.FrameRate = 10;
open(video2);

figure()
for i = 1: Size_mat_mesh_time
    title('App. error in diff. time')
    xlabel('Position')
    ylabel('App. error')
    ylim([0, max_error * 1.1])
    plot(Space_vect, abs(Mesh_mat(i,:) - sin(pi * Space_vect) * ...
        exp(- pi * pi * (i-1) * Time/ Time_step)))
    hold on
    pause(0.05)
    legend(strcat('Time',num2str((i-1)*Time/Time_step)))
    F2(i) = getframe(gcf);
%     fprintf('iteration %d',i)
end
saveas(gcf,'App_error_in_different_time','jpg')

writeVideo(video2,F2)
close(video2)

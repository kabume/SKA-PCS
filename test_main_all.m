% close all;
addpath data/
addpath data/img
addpath data/uvw
addpath data/vis

addpath lib/
addpath alg/

try
    run('../../../lib/irt/setup.m');
end

%% run parameters
% 0 - loads new data from file based on the dataset number supplied
% 1 - generates new data
% 2 - uses the data in matlab's workspace
gen_data = 1;
gen_figures = 1;
gen_only_average_figures = 0;
free_memory = 0;

save_dataset_number = 0; % number of the dataset to write files to
save_dataset_subnumber = 0; % number of the dataset to write files to

save_data_on_disk = 0; % flag
save_eps_files = 0; % flag
save_path = 'results/';
% save_path = '~/Data/aonose/results/';


num_tests = 1;
num_workers = 1; % number of tests to run in parallel; should be less than 
                 % the number of cores and is limited by the system memory for the variables


run_pdfb_bpcon_par_sim_rescaled = 1; % flag

run_pdfb_bpcon_par_sim_rescaled_precond = 1; % flag

run_pdfb_bpcon_par_sim_rescaled_precond_wave_par = 0; % flag

run_admm_bpconpar = 1; %flag

%% real data generation
use_real_visibilities = 1;
visibility_file_name = 'vla-3C129-pruned-extra.mat';

param_real_data.image_size_Nx = 512*2;
param_real_data.image_size_Ny = 512*2;
param_real_data.pixel_size = 1;

param_real_data.use_undersamplig = 0;

param_real_data.p = 10000;

param_real_data.fu_undersampling_type = 'uniform'; % 'uniform', 'gaussian', 'general-gaussian'

param_real_data.fu_g_sigma = pi/4; % variance of the gaussion over continous frequency

param_real_data.fu_ggd_beta = pi/4; % scale parameter for ggd
param_real_data.fu_ggd_rho = 0.07; % shape parameter for ggd

if use_real_visibilities % force only one test
    num_tests = 1;
    num_workers = 1;
end


%% simulated data generation
use_simulated_data = 0;
input_snr = 50; % noise level (on the measurements)
use_different_per_block_input_snr = 0;
per_block_input_snr_delta = 20; % delta defining the interval around the 
                                % input_snr from where the different snr levels are randomly chosen

image_file_name = 'CYGCBEST-1024.fits';

%% various config parameters
verbosity = 1;

nlevel = 4; % wavelet level
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

use_gridded_data = 0; % flag setting for generating gridded data

compute_block_op_norm = 0; % flag to compute the operator norm for each block

use_symmetric_fourier_sampling = 0;

%% definition for the stoping criterion
% options: 
% l2_ball_definition -> 'sigma', 'chi-percentile', 'value'
% stopping_criterion -> 'sigma', 'chi-percentile', 'l2-ball-percentage', 'value'

l2_ball_definition = 'value';
stopping_criterion = 'l2-ball-percentage';

param_l2_ball.stop_eps_v = [sqrt(2*307780)];
param_l2_ball.val_eps_v = 1.04 * param_l2_ball.stop_eps_v;

param_l2_ball.sigma_ball = 2;
param_l2_ball.sigma_stop = 2;

param_l2_ball.chi_percentile_ball = 0.99;
param_l2_ball.chi_percentile_stop = 0.999;

param_l2_ball.l2_ball_percentage_stop = 1.0001;

use_same_stop_criterion = 1; % forces the distributed criterion to be scaled
                             % such that same norm is imposed as in the nondistributed setup
 
%% sparsity prior
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used

         
%% samplig pattern parameters
% options 'gaussian', 'file', 'gaussian+large-holes', 'general-gaussian', 'gaussian-x2', 'file+undersample', 'gaussian+missing-pixels'

% sampling_pattern = 'file+undersample';
sampling_pattern = 'gaussian';
% sampling_pattern = 'general-gaussian';
% sampling_pattern = 'gaussian-x2';
% sampling_pattern = 'gaussian+large-holes';
% sampling_pattern = 'gaussian+missing-pixels';
% sampling_pattern = 'file';


param_sampling.p = 0.5; % number of measurements as proportion of number of pixels to recover

% file
% param_sampling.f_file_name = 'test_ska_240s.uvw.mat';
param_sampling.f_file_name = 'test_vla_10s.uvw.mat';



% 'file+undersample'
param_sampling.fu_file_name = 'test_vla_10s.uvw.mat';
param_sampling.fu_undersampling_type = 'gaussian'; % 'uniform', 'gaussian', 'general-gaussian'
param_sampling.fu_g_sigma = pi/4; % variance of the gaussion over continous frequency
param_sampling.fu_ggd_beta = pi/4; % scale parameter for ggd
param_sampling.fu_ggd_rho = 0.07; % shape parameter for ggd

% 'gaussian'
param_sampling.g_sigma = pi/4; % variance of the gaussion over continous frequency

% 'gaussian+missing-pixels'
param_sampling.gmp_hole_prob = 0.1; % probability of single pixel hole for 'gaussian+missing-pixels'
param_sampling.gmp_sigma = pi/4; % variance of the gaussion over continous frequency

% 'gaussian+large-holes'
param_sampling.gh_sigma = pi/4; % variance of the gaussion over continous frequency
param_sampling.gh_sigma_holes = pi/3; % variance of the gaussion for the holes
param_sampling.gh_hole_number = 8000; % number of holes to introduce for 'gaussian+large-holes'
param_sampling.gh_hole_size = pi/60; % size of the missing frequency data

% 'general-gaussian'
param_sampling.ggd_beta = pi/2; % scale parameter for ggd
param_sampling.ggd_rho = 0.3; % shape parameter for ggd

% 'gaussian-x2'
param_sampling.gx2_sigma_1 = pi/16; % variance of the gaussion over continous frequency
param_sampling.gx2_sigma_2 = pi; % variance of the gaussion over continous frequency
param_sampling.gx2_ratio = 0.995; % ratio of low freq data


%% block structure

regenerate_block_structure = 1;


param_block_structure.use_density_partitioning = 0;
param_block_structure.density_partitioning_no = 1;

param_block_structure.use_uniform_partitioning = 0;
param_block_structure.uniform_partitioning_no = 4;

param_block_structure.use_equal_partitioning = 1;
param_block_structure.equal_partitioning_no = 1;

param_block_structure.use_manual_frequency_partitioning = 0;
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, ...
    icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)

param_block_structure.use_manual_partitioning = 0;
param_block_structure.partition = [1000 2000 4000];

%% preconditioning

param_precond.gen_uniform_weight_matrix = 1; %set weighting type
param_precond.uniform_weight_sub_pixels = 1;

%% get input data
script_get_input_data;


%% PDFB parameter structure sent to the algorithm
param_pdfb.im = im; % original image, used to compute the SNR
param_pdfb.verbose = verbosity; % print log or not
param_pdfb.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb.nu2 = evl; % bound on the norm of the operator A*G
param_pdfb.gamma = 1e-6; % convergence parameter L1 (soft th parameter)
param_pdfb.tau = 0.49; % forward descent step size
param_pdfb.rel_obj = 0;1e-3; % stopping criterion
param_pdfb.max_iter = 20; % max number of iterations
param_pdfb.lambda0 = 1; % relaxation step for primal update
param_pdfb.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb.sol_steps = [inf]; % saves images at the given iterations

param_pdfb.use_reweight_steps = 0;
param_pdfb.use_reweight_eps = 0;
param_pdfb.reweight_steps = [200 inf];
param_pdfb.reweight_rel_obj = 1e-5; % criterion for performing reweighting
param_pdfb.reweight_min_steps_rel_obj = 100;
param_pdfb.reweight_alpha = 0.01;
param_pdfb.reweight_alpha_ff = 0.5;


%% PDFB parameter structure sent to the algorithm
param_pdfb_precond.im = im; % original image, used to compute the SNR
param_pdfb_precond.verbose = verbosity; % print log or not
param_pdfb_precond.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb_precond.nu2 = evl_precond; % bound on the norm of the operator A*G
param_pdfb_precond.gamma = 1e-6; % convergence parameter L1 (soft th parameter)
param_pdfb_precond.tau = 0.49; % forward descent step size
param_pdfb_precond.rel_obj = 0;1e-3; % stopping criterion
param_pdfb_precond.max_iter = 20; % max number of iterations
param_pdfb_precond.lambda0 = 1; % relaxation step for primal update
param_pdfb_precond.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb_precond.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb_precond.sol_steps = [inf]; % saves images at the given iterations

param_pdfb_precond.use_proj_elipse_fb = 1;
param_pdfb_precond.elipse_proj_max_iter = 200;
param_pdfb_precond.elipse_proj_min_iter = 1;
param_pdfb_precond.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

param_pdfb_precond.use_reweight_steps = 0;
param_pdfb_precond.use_reweight_eps = 0;
param_pdfb_precond.reweight_steps = [500:100:2000 inf];
param_pdfb_precond.reweight_rel_obj = 1e-5; % criterion for performing reweighting
param_pdfb_precond.reweight_min_steps_rel_obj = 50;
param_pdfb_precond.reweight_alpha = 0.01; % omega^(0)
param_pdfb_precond.reweight_alpha_ff = 0.5; % exponential decay factor for omega^(0)

%% ADMM parameter structure sent to the algorithm
param_admm.im = im; % original image, used to compute the SNR
param_admm.verbose = verbosity; % print log or not
param_admm.gamma = 1e-6*evl; % convergence parametter
param_admm.nu = evl; % bound on the norm of the operator A
param_admm.rel_obj = 0;1e-4; % stopping criterion
param_admm.max_iter = 20; % max number of iterations
param_admm.tight_L1 = 0; % indicate if Psit is a tight frame (1) or not (0)
param_admm.max_iter_L1 = 100; % max number of iterations to be performed for estimating the L1 prox
param_admm.rel_obj_L1 = 1e-3; % stopping criterion for the L1 prox
param_admm.pos_L1 = 1; % constrain for the positivity
param_admm.nu_L1 = 1; % bound on the norm of the operator Psi
param_admm.verbose_L1 = 0;  % print log or not
param_admm.sol_steps = inf; % saves images at the given iterations


%% compute the solution
fprintf('Starting algorithms:\n\n');
tstart = tic;

script_run_all_tests_serial;

tend = toc(tstart);
fprintf('All algorithms runtime: %ds\n\n', ceil(tend));

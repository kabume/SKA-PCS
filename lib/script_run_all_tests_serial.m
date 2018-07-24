%%

if run_pdfb_bpcon_par_sim_rescaled
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, ~, ~, result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, W, Psi, Psit, Psiw, Psitw, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%%
if run_pdfb_bpcon_par_sim_rescaled_precond
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
    result_st.xcorr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ~, ~, ...
            result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_precond runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond';
    param_structure_name = 'param_pdfb_precond';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


%%
if run_pdfb_bpcon_par_sim_rescaled_precond_wave_par
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
    result_st.xcorr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond_wave_par\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ~, ~, ...
            result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond_wave_par(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_precond_wave_par runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond_wave_par';
    param_structure_name = 'param_pdfb_precond';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


%%
if run_admm_bpconpar
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running admm_bpconpar\n');
        [result_st.sol{i}, ~, result_st.L1_v{i}, result_st.L2_v{i}, ...
            result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
            = admm_bpconpar(yT{i}, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), A, At, T, W, Psiw, Psitw, param_admm);
        tend = toc(tstart_a);
        fprintf(' admm_bpconpar runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'admm_bpconpar';
    param_structure_name = 'param_admm';

    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end





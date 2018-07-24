# SKA-PCS
It's the SKA summer internship project, astronomical image reconstruction using parallel compress sensing technique

error: op_nufft (line 25)
fatrix2 has some very weird bugs with subindexing; force st.p to be a (sparse) matrix

error: op_p_nufft (line 101)
            [~, ~, Gb, ~] = op_nufft([p{q, 1} p{q, 2}], N, Nn, No, Ns);

error: script_get_input_data (line 189)
        [A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);

error: test_main_all (line 196)
script_get_input_data;
note the line 294 in nufft_init.m, run with normal.
lib:
Matlab wavelet toolbox
Matlab Statistics and Machine Learning Toolbox Apps
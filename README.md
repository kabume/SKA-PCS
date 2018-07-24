# SKA-PCS
It's the SKA summer internship project, astronomical image reconstruction using parallel compress sensing technique
错误使用 op_nufft (line 25)
fatrix2 has some very weird bugs with subindexing; force st.p to be a (sparse) matrix

出错 op_p_nufft (line 101)
            [~, ~, Gb, ~] = op_nufft([p{q, 1} p{q, 2}], N, Nn, No, Ns);

出错 script_get_input_data (line 189)
        [A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);

出错 test_main_all (line 196)
script_get_input_data;
注释掉294行后运行正常

wavelet toolbox
Statistics and Machine Learning Toolbox Apps
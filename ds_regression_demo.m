% Set bins for regressors corresponding to each widefield frame
wf_regressor_bins = [wf_t;wf_t(end)+1/wf_framerate];

% Create regressors
stim_regressor = histcounts(stimOn_times,wf_regressor_bins);
move_regressor = histcounts(stim_move_time,wf_regressor_bins);

regressors = {stim_regressor;move_regressor};

% Set time shifts for regressors
t_shifts = {[-5:30];[-30:30]};

% Set cross validation (not necessary if just looking at kernels)
cvfold = 5;

% Do regression
[kernels,predicted_signals,explained_var,predicted_signals_reduced] = ...
    ap.regresskernel(regressors,wf_V,t_shifts,[],[],cvfold);

% Convert kernels V to pixels
kernels_px = cellfun(@(x) plab.wf.svd2px(wf_U,permute(x,[3,2,1])),kernels,'uni',false);

% Plot kernels
plot_kernel = 1;
ap.imscroll(kernels_px{plot_kernel});
axis image;ap.wf_draw('ccf');
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG',[],1.5));

plot_kernel = 2;
ap.imscroll(kernels_px{plot_kernel});
axis image;ap.wf_draw('ccf');
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG',[],1.5));










%% Regression demo
%
% This demos how to do linear regression with the ap.regresskernel
% function. (There are other ways to do this in MATLAB that I'm not familar
% with, so it may be worth exploring other built-in functions to do this).
%
% This is useful when trying to disentangle responses to multiple
% overlapping events (like stimuli and movements), or when trying to
% determine the function linking one pattern of activity to another (e.g.
% two simultaneously recorded regions).
%
% Linear regression determines a kernel that, when convolved with the
% regressors, produces a predicted signal with a least-squares error to
% your original signal.
%
% In this example, we will try to disentangle stimulus responses and
% movement responses, even though both happen close together and are mixed
% in the average stimulus-aligned response:
%
% - our signal is the widefield data (this is what we want to predict).
% - our regressors will be binary vectors of where stimuli and movements
% start (this is what we will use to predict our signal).
% - the kernels we will find are the de-mixed responses to stimuli and movements.


%% Example dataset

% Load this example widefield dataset: 
animal = 'AP015';
rec_day = '2024-01-25';
rec_time = '1650';
verbose = true;

load_parts.widefield = true;
ap.load_recording;

% This is from a mouse in early learning of the sensorimotor task. 
% The relevant variables from this dataset are: 
%
% - wf_V/wf_U/wf_t = widefield V/U/timepoints
% - stimOn_times = onset of stimuli
% - stim_move_time = start times of movements after stimuli
% - stim_to_move = reaction time from stimulus to movement time

%% Mixed response

% The question we want to ask in this data is: what does the response to
% the stimulus look like? 

% Stimuli and movements occur very close together in this data, so
% averaging our data to the stimulus onset will include both stimulus- and
% movement-related responses.

% [EXERCISE]
% 1) Make an average widefield movie aligned to the stimulus onset. Can you
% tell what part of activity is related to stimulus, and what's related to
% movement?
%
% 2) Draw an ROI over the visual cortex, then plot the activity of that ROI
% aligned to the stimulus for each trial, sorting trials by reaction time.
% So: y-axis is trial sorted by reaction time, x-axis is time relative to
% stim onset. Can you see where the stim- and movement-related activity is?

%% Regressing stimulus and movement responses

% You can do linear regression using `ap.regresskernel`. Look through the
% documentation (help ap.regresskernel) for the syntax.

% Because linear regression is a linear operation, and all linear
% operations can be performed either in pixel space or V space, we can
% perform regression directly on the V's and save ourselves lots of space
% and time. 
%
% So - the "signals" that we are trying to predict are the V's. This
% dataset has saved 2000 components, but we can get a pretty accurate
% prediction using only a fraction of these (e.g. 200) to save us further
% time and space. 
%
% Let's define the number of components to use as 200: 

use_n_components = 200;

% The first step is to set up your regressors - this is what will be used
% to predict your data. In this example, we will make two regressors:
% stimulus onsets, and movement onsets (post-stimulus movements only).
%
% The time dimension of the regressors has to match the data, because it is
% a one-to-one mapping. For example, if your data has 5 timepoints, each
% regressor also needs to have 5 timepoints. 
%
% So for each regressor, we want a binary vector which is the same length
% as the number of widefield frames, with 0's when an event did not happen
% and a 1 when it did. For example, if you have 10 frames, and the stimulus
% was delivered during the 4th frame, you want this as your regressor: 
% [0 0 0 1 0 0 0 0 0 0]

% [EXERCISE]
% Create type of regressor for the stimulus and movement onsets. Stim
% onsets are in `stimOn_times`, and post-stim movement onsets are in
% `stim_move_times`. Note these are timestamps of events, so you'll need to
% bin them to match the widefield frames. Do this using the `histcounts`
% function and `wf_t`, being the widefield frame timestamps.

stim_regressor = []; % make the stim regressor here
move_regressor = []; % make the move regressor here

% We now package these regressors into a cell array to input into the
% function later: 
regressors = {stim_regressor,move_regressor};

% Our regressors currently only correspond to the onset of each event, but
% we want kernels that span some timecourse, because the activity follows
% some pattern across time. 
%
% We can do this by shifting our current regressors forwards and backwards
% in time. So our current regressors give us the activity when lag = 0s,
% shifting them forward one timestep will give us activity when lag = +1
% frame, etc. We don't have to make these lagged regressors by hand, the
% function does that for you, by specifying the time lags to use.
%
% These time shifts can be specified for each regressor, and should
% correspond roughly to how long we expect the response to be. 
% Let's set the stimulus regressor to have lags of 0-30 frames (assumes the
% stimulus response is ~1s and starts from stimulus onset), and the
% movement regressor to have lags of -30:30 frames (allows the movement
% response to be longer, and allows movement-related activity to start
% before the actual movement starts):
time_shifts = {0:30,-30:30};

% We can now do the regression with `ap.regresskernel`: 
[kernels,predicted_signals,explained_var,predicted_signals_reduced] = ...
    ap.regresskernel(regressors,wf_V(1:use_n_components,:),t_shifts);


%% Working with regression outputs

% kernels: 
% The kernels, convolved with the regressors and summed
% together across regressors, best estimate the data. In this case, they
% correspond to the "stimulus response" and the "movement response".

% [EXERCISE]
% Our kernels are currently given as V's. Convert these to pixels, and plot
% the movies. 
% Do these make sense? 
% How do they compare to the average stimulus-aligned movie you made above?

% predicted_signals (_reduced)
% These are the signals predicted by convolving the kernels with the
% regressors. `predicted_signals` are using all regressors,
% `predicted_signals_reduced` are using all regressors *except* for a given
% regressor. So while `predicted_signals` has size n_components x time,
% `predicted_signals_reduced` has size n_components x time x n_regressors,
% where in this case we have 2 regressors (stim and movement). The first
% page of this (`predicted_signals_reduced(:,:,1)`) is the predicted signal
% *without* the first regressor, in other words: this is the predicted
% signal using only movement. The same applies to the second page and
% movement. 
%
% We can use these predicted signals to give us "residuals", which is the
% activity left over from regressing out a particular component. For
% example, if we want the activity that is just related to movement, we
% can subtract the activity predicted from all kernels except for the
% movement kernel. In our case, that's the 2nd kernel, so the second page
% of the `predicted_signals_reduced`: 
move_activity = wf_V(1:n_components,:) - predicted_signals_reduced(:,:,2);

% [EXERCISE]
% Find activity which is just related to the stimulus by subtracting
% activiy predicted by all kernels besides the stimulus. Draw an ROI over
% the visual cortex, and plot this residual activity as a heatmap as above:
% y-axis is trial sorted by reaction time, x-axis is time relative to
% stimulus onset. How does this compare to the heatmap you made earlier
% from the raw data? Does that make sense? 

%% Linear regression caveats

% Linear regression makes an assumption that is almost never true: that a
% given type of response is always the same. This is particularly not true
% for movement, e.g. the more vigorous the movement, the stronger the
% activity. What regression does in this case is find the best-fit
% "movement response" across all movement responses. BUT, because it's not
% a one-kernel-fits-all situation, it will not always do a great job of
% "removing" movement-related activity. 

% Also, the less jitter there are between events, the harder it will be to
% separate them. For example, if there is always exactly 0.5s between two
% events, it will be impossible to distinguish where each particular
% response starts. 

% [EXERCISE]
% Try loading this data, and performing the same type of regression above: 
animal = 'AP015';
rec_day = '2024-02-01';
rec_time = '1717';
verbose = true;
%
% Make a movie of the stimulus and movement kernels. How are they different
% from the kernels you made in the previous example? Why is it different?

%% (to be included later)

% Not included in this demo yet are these other features: 
% - when and how to include a constant as a regressor
% - ridge regression using a penalty value lambda
% - cross-validation and explained variance to examine goodness-of-fit







%% Electrophysiology demo 1
%
% A demo on how to work with electrophysiology data
% (do plab_data_demo first)
% Whenever you see a line of code, run it and continue reading
% [EXERCISE] marks a section with instructions to complete

%% Github repositories to download

% Clone these repos and add to the Matlab path 
% (Home > Set Path > Add with Subfolders...)
%
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab
% https://github.com/kwikteam/npy-matlab
% https://github.com/petersaj/neuropixels_trajectory_explorer

%% Neuropixels probes

% We use Neuropixels probes to record electrophysiological activity. 

% There are two types of Neuropixels probes: 
%
% Neuropixels 1.0 (one shank): 
% https://www.nature.com/articles/nature24636
%
% Neuropixels 2.0 (4-shank): 
% https://www.science.org/doi/full/10.1126/science.abf4588

% Both types of Neuropixels probes are 1 cm long, and include 1000 sites
% along the length of the shank. We can record from 384 channels at the
% same time (but one of these channels is a "reference" channel that
% doesn't record voltage, so we collect data from 383 channels). 

% We can also choose which sites to record from. In a 1-shank probe, this
% is almost always the 384 channels closest to the tip. In a 4-shank probe,
% this is more configurable, e.g. the bottom 96 channels of all 4 shanks,
% 192 channels from 2 shanks, etc. 

% A pretty comprehensive list of Neuropixels resources and tools is here:
% https://github.com/Julie-Fabre/awesome_neuropixels

% If you want to watch demos about Neuropixels, a good place to start is
% the UCL neuropixels course, which places all recorded lectures online.
% You can find the latest course here: 
% https://www.ucl.ac.uk/neuropixels/courses

%% Preprocessing

% There's no need to download this code, but for information:

% We use Kilosort 4 to detect and sort spikes: 
% https://github.com/MouseLand/Kilosort

% Kilosort can detect artifacts that are not action potentials, so we
% identify and remove artifacts using Bombcell:
% https://github.com/Julie-Fabre/bombcell


%% Example dataset

% Run these lines to load example data for the following sections
% This data contains units from the visual cortex (as well as underlying
% hippocampus).

animal = 'AP003';
rec_day = '2023-06-07';
rec_time = '1542';

animal = 'AP019';
rec_day = '2024-05-07';


verbose = true;
ap.load_recording;


%% Templates

% Kilosort works by template matching, which looks for repeated patterns of
% activity across the probe that correspond to the spike from a given
% neuron. These templates are stored in the variable 'templates', with size
% templates x sample x channels
% 
% The sample rate is 30kHz, so each sample = 1/30 ms
% The "spike time" corresponds to sample 21 in the template.
%
% Let's look at an example template. This will show the waveform on each
% channel, which constitutes the entire template. In other words: whenever
% this neuron spikes, this is what the spike looks like across the probe:
example_template_num = 169;
example_template = permute(templates(example_template_num,:,:),[3,2,1]);
figure;imagesc(example_template);
xlabel('Sample');
ylabel('Channel');
title(sprintf('Template %d',example_template_num));

% We can also plot the waveform of each channel overlaid to better see
% the shapes of the waveforms:
figure;plot(example_template');
xlabel('Sample');
ylabel('~mV (not exactly)')
title(sprintf('Template %d',example_template_num));

% The words "template", "unit", and "neuron" can be used interchangably,
% depending on where they're used. Often, "unit" is used instead of
% "neuron", because we can't be 100% sure that we're perfectly separating
% the spikes from one neuron (e.g. one "unit" might have a few spikes from
% a different neuron). If we think our unit contains multiple neurons, or
% we group multiple units together, we call that "multiunit" activity, as
% opposed to "single unit", which we think are single neurons. 

% [EXERCISE] 
% We typically define "the waveform" of a unit as the waveform on the
% channel with the largest amplitude (which could be positive or negative).
% Find and plot "the waveform" for the example template above. What channel
% does that correspond to?
%
% The load script ap.load_ephys calculates this already as 'waveforms',
% which is size template x sample. Check that your waveform matches the one
% from this variable. 

% The position of each channel (aka 'site') is stored in
% 'channel_positions'. A map of the channels in space can be seen by
% plotting these positions:
figure;
plot(channel_positions(:,1),channel_positions(:,2),'.k')
set(gca,'YDir','reverse');
xlabel('Probe X (\mum)');
ylabel('Probe Y (\mum)')
xlim([0,50]);
% In this variable, 0 is the top, and 3840 is near the tip (the deepest
% point). Try zooming into this plot to see the arrangement: note that
% there are two columns of sites, so that there are 2 recording positions
% at each depth.

% Where is each unit along the probe? Since templates span the entire
% probe, we define this as a weighted average of the waveform amplitude for
% the depth position of each channel (e.g. a big waveform on the channel
% 100 depth=2840 and a small waveform on the channel 101 depth=2820 could
% give a unit depth = 2835).
%
% The load script ap.load_ephys calculates this as 'template_depths'. Here
% is a plot of the depth for each unit:
figure;plot(template_depths,'.k');
xlabel('Template');
ylabel('Depth (/mum)');

% [EXERCISE] 
% Plot a histogram of the template depths in 100um bins. Note the pattern:
% depending on where the probe goes through, it picks up different numbers
% of units in different areas. Also note: the probe on this day was only
% ~2600um into the brain, which means there should not be any neurons with
% template depth less than 1240um (3840um full length - 2600um in brain =
% 1240um outside the brain), but we see there are some. This means our
% quality control is not perfect - we are detecting "spikes" that are not
% actually neurons.

% In summary, here are the relevant template variables: 
% - templates
% - template_depths
% - waveforms


%% Spikes

% Whenever Kilosort finds an instance of a template in the raw data, it
% logs that sample as a "spike" for that template/unit/neuron. During
% preprocessing, we convert that sample into a time, which is loaded in as
% 'spike_times_openephys' (Open Ephys is the recording program we use).
% In order to synchronize these times with things we care about, the load
% script ap.load_ephys converts times on the recording clock into Timelite
% times. (Just for reference - it does this using the "flipper" signal,
% which is a randomly flipping signal that goes into both our
% electrophysiology recording and Timelite). 
%
% The resulting vector is 'spike_times_timelite', which is size spikes x 1.
% In other words, this represents the time (on the Timelite clock) for
% every spike detected from every template across the entire recording.

% We can identify which template each spike belongs to with the variabe
% 'spike_templates'. This variable is also size spikes x 1, so each spike
% time has a corresponding template identity. 
%
% For example, the nth spike detected happens at 'spike_times_timelite(n)',
% and it belongs to the template 'spike_templates(n)':
use_spike = 1000;
fprintf('Spike %d happened at time %gs and belonged to template %d\n', ...
    use_spike,spike_times_timelite(use_spike),spike_templates(use_spike));

% As an example for how to use these variables, let's say we want to know
% how many spikes were detected for template 10: 
use_template = 10;
use_template_spike_n = sum(spike_templates == 10);
fprintf('Template %d had %d spikes\n',use_template,use_template_spike_n);
% 

% Let's say we want to plot how many spikes occurred for this unit across
% the recording: 
use_spike_times = spike_times_timelite(spike_templates == use_template);
figure; histogram(use_spike_times);
xlabel('Recording time (s)');
ylabel('Number of spikes per bin');
title(sprintf('Unit %d',use_template));
xline(0,'r');

% Look at the time (x) axis in the histogram plot above, and note that is
% starts at a negative time, and zero is somewhere in the second half
% (where the red line is drawn). This is because the spike times are
% relative to the Timelite for this particular recording, where zero is the
% start of this recording. Since this recording was not the first one of
% the day (there were two others before it), the negative times correspond
% to spikes that happened before this particular recording.

% [EXERCISE] 
% Calculate the number of spikes for each template in the whole day. Plot
% the number of spikes against the depth of each template, to make a plot
% that has template depth on the y-axis and number of spikes on the y-axis.
%
% Next, make the same plot, but calculating the number of spikes only in
% the loaded recording (hint: spike times should be within Timelite times).
% 
% Bonus points if you calculate spike number using the 'accumarray'
% function. This function can be confusing to start working with, but can
% be powerful once you know how to use it. 

% We saw in the templates section above that ap.load_ephys calculates the
% depth for each template. This is also stored on a spike basis by
% ap.load_ephys in the variable 'spike_depths', which has size spikes x 1
% like the other spike variables. In other words: the nth spike was
% positioned at 'spike_depths(n)':
use_spike = 1000;
fprintf('Spike %d was positioned at depth %gum\n', ...
    use_spike,spike_depths(use_spike));

% [EXERCISE] 
% Calculate and plot the number of spikes in 500um segments along the probe
% (e.g., how many spikes were detected on the probe between 1500-2000um?)

% In summary, here are the relevant spike variables: 
% - spike_times_timelite
% - spike_templates
% - spike_depths


%% Binning spikes

% It's usually easier to work with spikes when they're grouped into time
% bins, rather than working with all spike times independently. 

% Matlab's main function to bin data is 'histcounts'. Here's an example
% use, where we bin the spikes from an example unit into 10s time bins: 
use_template = 10;
use_spikes = spike_times_timelite(spike_templates == use_template);

time_bin_size = 10;
time_bins = min(use_spikes):time_bin_size:max(use_spikes);

use_spikes_binned = histcounts(use_spikes,time_bins);

% The x-value for each bin is the center of each bin. So if the bin edges
% are represented by |, and the centers by *, that would look like: 
% |  *  |  *  |  *  |
% Note that for n bins (3 *'s), there are n+1 bin edges (4 |'s). To
% calculate bin centers, we take (1:end-1) bins + diff(bins)/2 (halfway
% between each bin): 
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% We usually care about the rate (spikes per second), rather than the
% absolute number of spikes. We calculate that by dividing the number of
% spikes in each bin by the bin size:
use_spikes_binned_rate = use_spikes_binned/time_bin_size;

% Finally, lets plot the spike rate over time:
figure;plot(time_bin_centers,use_spikes_binned_rate);
xlabel('Time (s)');
ylabel('Spikes/s');

% [EXERCISE] 
% Just as we binned by time above, we can also bin by space along the
% probe.
% 
% Make a 2D variable, where bins in dimension 1 are 50um increments on the
% probe, and bins in dimension 2 are 500ms increments.
% 
% Plot that matrix, and underneath that axis, plot the wheel velocity
% ('wheel_velocity'). Link the two axes in 'x' with 'linkaxes', so that
% when you zoom into one plot, it zooms into the other. 
% 
% What's the relationship between spikes along the probe and wheel
% movement?

% [EXERCISE] 
% The variable you made above with spikes binned by time and space gives
% the multiunit activity (spikes from many neurons) at different depths
% along the probe. Try using the function 'corr' on this matrix to find the
% correlation of each depth with all other depths (note: be careful about
% the orentation of your input to 'corr'). What do you notice? 


%% Peri-stimulus time histograms (PSTH)

% The most common way to analyze spiking activity is by aligning spikes to
% events. This is done by binning spikes in time bins as above, but
% relative to an event rather than across the whole recording. 

% The loaded recording contains passive presentations of visual stimuli on
% the left, center and right screen. Our probe is in the left visual
% cortex, so we expect to see responses to stimuli on the right screen. 
% 
% Let's calculate the PSTH for an example neuron to stimuli in the right
% screen: 

% Using Bonsai events, we get the X (azimuth) stimulus positions:
stim_x = vertcat(trial_events.values.TrialStimX);

% Using Timelite, find the stim times when x = +90 (right screen)
right_stim_times = stimOn_times(stim_x == 90);

% Create time bins with a certain width around each stim
% (define bin size in seconds)
bin_size = 0.001; 
% (define the start and end of the window, in seconds)
psth_window = [-0.5,1];
% (define the relative bins for the PSTH, edges and centers (to plot))
psth_bins = psth_window(1):bin_size:psth_window(2); 
psth_bin_centers = psth_bins(1:end-1)+ diff(psth_bins)/2;
% (get the bins around every stim: stim is on dim 1, bin is on dim 2)
stim_bins = right_stim_times + psth_bins; 

% (select a unit, loop through each stimulus, bin the spikes)
use_unit = 180;
use_unit_spikes = spike_times_timelite(spike_templates == use_unit);

unit_psth = nan(length(right_stim_times),length(psth_bins)-1);
for curr_trial = 1:length(right_stim_times)
    % (get the binned spikes for this trial, divide by the bin size to get
    % spike rate)
    unit_psth(curr_trial,:) = ...
        histcounts(use_unit_spikes,stim_bins(curr_trial,:))/bin_size;
end


% all_counts = cell2mat(arrayfun(@(x) histcounts(use_unit_spikes, stim_bins(x,:)), (1:length(right_stim_times))', 'UniformOutput', false));


% The bin size is currently very small: 1ms, so small that each bin has
% either 0 or 1 spikes. If we plot the PSTH now, it will have high temporal
% precision, but it will look messy: 
figure;
plot(psth_bin_centers,mean(unit_psth));
xlabel('Time from stim');
ylabel('Spikes / s')
xline(0,'r');
title(sprintf('Unit %d: PSTH to right stim',use_unit));

% To clean this up, one thing we can do is apply a gaussian filter to
% smooth our data, which loses some temporal precision but provides a more
% easily interpretable plot. There are a few ways to apply a gaussian
% filter, this way uses the function 'smoothdata', and shows the effect of
% smoothing with different sized windows:
figure; hold on
smooth_windows = 1:20:100;
for curr_window = 1:length(smooth_windows)
    plot(psth_bin_centers,smoothdata(mean(unit_psth), ...
        'gaussian',smooth_windows(curr_window)));
end
xlabel('Time from stim');
ylabel('Spikes / s')
xline(0,'linewidth',2)
title(sprintf('Unit %d: PSTH to right stim',use_unit));
legend(cellfun(@(x) sprintf('Smoothing window: %d',x), ...
    num2cell(smooth_windows),'uni',false));

% Try zooming into the plot and noticing the differences between the
% different smoothing windows, particularly the trade-off between
% smoothness and temporal precision. Which value do you think works best?

% [EXERCISE] 
% In the above example, we had very small bins (1 ms), which we then
% smoothed. How does this compare to just using large bins and not
% smoothing? In other words, what's the difference between 1 ms bins
% smoothed with a 30ms window, vs. just using 30ms bins? Try comparing
% PSTHs with this difference using 10,20, and 50 ms windows

% [EXERCISE] 
% The above example used a for loop to bin spikes for each trial. Try
% binning spikes for each stimulus using the function 'arrayfun', which can
% act as a one-line for loop. 


%% Determining responsive cells

% We often want to extract cells that do something we're interested in. For
% example, in this case, we may want to pull out cells that respond to the
% visual stimulus for further analysis. 

% There are many ways to do this, which vary in their sensitivity and
% features they emphasize. One method is to compare a "baseline" time to an
% "event" time, and choose cells which are significantly more active during
% the event compared to baseline.

% [EXERCISE] 
% Find units that are responsive to right-side stimuli in these steps: 
%
% 1) For each unit, count the number of spikes (histogram with one bin)
% 150ms before each stimulus (baseline) and 150ms after each stimulus
% (event). Store these counts in a matrix which is size units x 2
% [baseline,event] x stimulus.
%
% 2) For each unit, do a Wilcoxon signed rank test ('signrank' function) to
% get a p value comparing the baseline spikes to the event spikes.
%
% 3) Define "responsive" units as having p < 0.05 from the test above
%
% 4) Create an average PSTH around the stimulus for all units
%
% 5) Create a heatmap of PSTHs for responsive units, and another for
% non-responsive units
%
% 6) The above PSTHs may have a lot of variability due to baseline firing
% rates. We don't care much about baseline firing rates here, instead, we
% care about relative firing rates (e.g. when the stim came on, the spike
% rate doubled). Try normalizing the PSTHs by defining a baseline firing
% rate for each unit (maybe averaging 200ms before stim onset). Then
% calculate the normalized firing rate by subtracting and dividing by the
% baseline as normalized_rate = (rate-baseline_rate)/baseline_rate. Plot
% the heatmap of these noramlized PSTHs - what's the difference with the
% non-normalized version? What are the units for this plot (e.g. what does
% "2" mean)?


%% PSTH viewer

% I have a tool for viewing PSTHs: 'ap.cellraster'. This function pulls
% certain variables from the workspace (e.g. spike times, templates,
% depth), so data first has to be loaded in using the conventions in
% 'ap.load_ephys' (which is called by 'ap.load_experiment'). In other words
% - if you don't load data with 'ap.load_experiment', you may get an error
% about missing variables. 

% 'ap.cellraster' takes two inputs: 
% align_times - times of events to align (this can be a vector, e.g.
% stimulus times, or a cell array, e.g. cell 1 = stimulus times, cell 2 =
% reward times)
%
% align_groups - categories of entries in 'align_times' (e.g. if there are
% two alternating stimuli displayed 0.5s apart, the 'align_times' may be
% [0,0.5,1,1.5] and the 'align_groups' would be [1,2,1,2]).

% Here is an example usage to align activity to stimuli, using stimulus
% onset times as the alignment times, and stimulus azimuth as the grouping:
stim_x = vertcat(trial_events.values.TrialStimX);

stim_x = vertcat(trial_events.values.StimFrequence);

ap.cellraster(stimOn_times,stim_x);

% Here's a desription of the interface: 
% - left plot: a dot for each template, with depth
% on the y-axis and normalized rate on the x-axis (highest firing rate
% cells are on the right, lowest are on the left). 
%
% - center plot: waveforms of the template on a selection of channels
% (channels with the largest amplitude)
%
% - right plot: top is the PSTH, bottom is a raster plot (one dot for each
% spike, time on the x-axis, trial on the y-axis)
%
% - bottom plot: template amplitude across time (this is used to see if a
% unit drifts in or out of the recording over time, which can happen with
% subtle movement of the probe/brain). 

% [EXERCISE] 
% Try using the main functions of cellraster: 
% - click on dots in the left plot - this selects a unit to plot
% 
% - press space bar: this toggles between plotting all events, vs.
% grouping events by the grouping variable (in this case, stim azimuth).
% The default colors are blue/black/red if the grouping variables have both
% negative and positive numbers. In this case, the colors correspond to the
% stimulus positions blue = -90 (left), black = 0 (center), red = 90
% (right). What's the difference between these groups? Why?
%
% - press the up/down keys: this cycles through templates by depth
%
% - press 'm' (for 'multiunit'), and click on the left plot around depth
% 1500 and then again around depth 2500: this feature groups the spikes of
% all units within the selected depth. Since there are too many spikes to
% display in a coherent raster, the raster plot turns into a heatmap.
%
% - press 'u' (for 'unit'), and enter 152: This goes to unit 152. 
%
% press 't' (for 'time'), and enter [-0.1,0.5]: this changes the displayed
% time of the PSTH and raster (the default when you open is [-0.5,2]). 


% [EXERCISE] 
% Try calling 'ap.cellraster' with a second grouping variable. Even though
% the loaded recording was a passive protocol, the mouse did move the wheel
% occasionally.
%
% First, get the times at the start of each movement (use the
% 'wheel_move' variable). 
%
% Second, find the duration for each movement (offset time - onset time).
% Get the sort index of movement durations (e.g. the shortest movement is
% 1, the next shortest movemnet is 2). We'll use this as our grouping
% variable (which doesn't have to be categorical groups - it can also be a
% sorting).
%
% Call 'ap.cellraster' with a cell array for align times and grouping
% variables: the first cell should have the stim times and groups used
% above, the second cell should have movement times and sort index.
%
% Use 'm' (multiunit mode) to select the visual cortex, which is depth
% ~1500-2500. Press right/left to toggle between alignments (note at the
% top of the right plots, it switches betwen "Align 1" and "Align 2", which
% is just the order of alignments put into the function). On Align 2
% (movement), press space bar, which will sort the trials by movement time
% (note the y-axis of the heatmap switches between "Trial" and "Sorted
% trial"). 
%
% This is visual cortex - what do you notice about responses during visual
% stimuli vs. movement?


%% Planning trajectories and estimating live probe location

% I built a tool called the Neuropixels Trajectory Explorer to plan
% Neuropixels trajectories, which also connects to the manipulators during
% an experiment to show estimated live probe location. This is used by many
% labs, so I regularly answer questions online and fix bugs / make updates.

% This uses the Allen Common Coordinate Framework mouse brain atlas
% (shortened to "CCF" for "common coordinate framework", or sometimes the
% "ARA" for "Allen reference atlas"), which is described in this paper:
% https://www.sciencedirect.com/science/article/pii/S0092867420304025

% The repository is here: 
% https://github.com/petersaj/neuropixels_trajectory_explorer

% You can watch a video demo here from the UCL Neuropixels course: 
% https://www.youtube.com/watch?v=54VHDqzowwY&ab_channel=MatteoCarandini

% Open the program with this command: 
neuropixels_trajectory_explorer

% Documentation for how to use it is on the Github wiki page: 
% https://github.com/petersaj/neuropixels_trajectory_explorer/wiki
% 
% Have a look through "General use" and "Example trajectory planning"

% [EXERCISE] 
% Following the "Example trajectory planning", plan a trajectory with a
% 1-shank probe (Neuropixels 1.0) between the visual cortex ("Primary
% visual area") and the lateral geniculate nucleus ("Dorsal part of the
% lateral geniculate complex").
%
% Now add a second, 4-shank (Neuropixels 2.0) probe, and position it to
% record simultanously from the striatum ("Caudoputamen"), GPe ("Gobus
% pallidus external segment"), and GPi ("Globus pallidus interal segment").
%
% Save this trajectory by pressing Save > Save positions.

% [EXERCISE] 
% Load the positions saved in the above exercise. the variable
% "probe_areas" is a N probes x 1 cell array, with each entry containing
% CCF information for areas along the length of the trajectory inside the
% brain. The area name is in the field "name", the location of that area on
% the probe (relative to 0 = insertion point) is given by the field
% "probe_depth", and the shank is given by "probe_shank".
%
% On the 1-shank probe: what depth does the visual cortex span? 
%
% On the 4-shank probe: which shanks record from the striatum, GPe, and
% GPi?

% For information, there is another type of tool like this called Pinpoint,
% built in Nick Steinmetz's lab and modeled after the NTE. That is accessed
% on a browser here:
% https://data.virtualbrainlab.org/Pinpoint/
%
% If you like, play around with that program, and see which one you like
% better / what the pros and cons are. Since that program is the focus of
% someone in their lab, it might replace the NTE some day.

%% Next demo

% To learn how to align electrophysilogical data to histology and determine
% where your record was in the brain, stay tuned for the next demo.






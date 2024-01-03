%% Data demo
%
% A demo on how to use the server and load data
% Whenever you see a line of code, run it and continue reading
% [EXERCISE] marks a section with instructions to complete

%% Github repositories to download

% Clone these repos and add to the Matlab path (Home > Set Path > Add with
% Subfolders...)
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab

%% File structure and finding recordings

% To connect to the server: 
% - Open up the file explorer, select This PC in the left hand menu 
% - Right click on this PC and select map network drive 
% - Enter the following address: \\qnap-ap001.dpag.ox.ac.uk\APlab\ 
% - It should prompt you for credentials, please enter your MSD details in the following format: 
% - Username – MSD\(MSD username) 
% - Password – normal MSD password 

% The server contains these top folders:
% \Data: for all raw and preprocessed data (all users can write)
% \Users: all users have a folder for their own use (given users can write)
% \Lab: files for the general lab and rigs (only Andy can write)

% Data is organized on the server in these folders: 
% \Data
% |- animal name (initials & number, e.g. AP012)
% | |- recording date (yyyy-mm-dd, e.g. 2023-03-16)
% |  |- recording time (Recording_HHMM, e.g. Recording_1021 started at 10:21am)
% |  | |- recording data by modality (e.g. mousecam, widefield)
% |  |- day data spanning multiple recordings (e.g. ephys, widefield)
% |- non-recording data (e.g. histology)


%% Finding data from recordings

% The function 'plab.find_recordings' is used to find data from recordings
% The syntax is: plab.find_recordings(animal,recording_day,workflow)
% (type 'help plab.find_recordings' for documentation).
% 'animal' must be filled in, other specifics can be left out to return all
% relevant recordings.

% --- Working with plab.find_recordings

% This line finds all recordings done with animal AP010:
% (day and workflow are not entered)
recordings = plab.find_recordings('AP010');

% 'recordings' is a structure with one entry for each day containing:
% - day: day of recording
% - index: index of recording within the day for that animal (e.g. 3rd recording = [3])
% - recording: time of recording (HHMM)
% - workflow: Bonsai workflow (i.e. stimuli or task protocol)
% - mousecam: if there was mouse camera recording (0/1 = no/yes)
% - widefield: if there was widefield recording (0/1 = no/yes)
% - ephys: if there was electrophysiology recording (0/1 = no/yes)

% [EXERCISE]
% 1) How many days were recordings performed on AP010?
% 2) How many recordings were performed on AP010 across all days?
% 3) How many unique Bonsai workflows were used, and what were they?

% This line finds all recordings of AP010 on 2023-08-24: 
% (workflow is not entered)
recordings = plab.find_recordings('AP010','2023-08-24');

% [EXERCISE]
% Which recording modalities were recorded this day?

% This line finds all recordings of AP010 with workflow 'lcr_passive':
% (date is not entered)
recordings = plab.find_recordings('AP010',[],'lcr_passive');

% [EXERCISE]
% 1) How many recordings of this workflow included widefield? 
% 2) How many recordings of this workflow included electrophysiology? 
% 3) Make a new variable which is a subset of 'recordings' with widefield

% This line finds recordings of AP010 with workflow 'lcr_passive' on
% 2023-08-16 (all fields are entered)
recordings = plab.find_recordings('AP010','2023-08-16','lcr_passive');

% [EXERCISE]
% This day has 2 recordings of the same workflow. This usually is
% because there was an issue with the first one and it needed to be re-run.
% In the number order of that day's recordings, which ones were
% 'lcr_passive'? (hint: 'index')

% The 'workflow' can include * as a wildcard. For example, in the task
% workflow 'stim_wheel_right', there is a 'stage1' and 'stage2' version.
% This returns specifically recordings with stage1:
recordings = plab.find_recordings('AP010',[],'stim_wheel_right_stage1');
% And this returns recordings with either stage: 
recordings = plab.find_recordings('AP010',[],'stim_wheel_right*');

% [EXERCISE]
% Write code to return the day AP010 switched from stage1 to stage2.


%% Constructing server paths

% The class 'plab.locations' is used to find or construct paths on the server.

% --- Standardized generic locations

% This is the general path to the server:
plab.locations.server_path
% This is the data path on the server:
plab.locations.server_data_path
% This is where local data should be kept on each computer: 
plab.locations.local_data_path

% [EXERCISE]
% Use the function 'fullfile' and the above information to generate the
% path: server/Users/(your name)/test_path

% --- Recording locations

% The method 'plab.locations.filename' is used to construct paths/filenames
% for recording data
% (type 'help plab.locations.filename' for documentation). 
% The syntax is:   constructed_filename = filename('server' or 'local',animal,rec_day,rec_time,folder1,...,folderN)

% For example: 
% This constructs the path for animal AP001 on day 2000-01-01 at 12:00:
plab.locations.filename('server','AP001','2000-01-01','1200')
% Input arguments can be added or removed as necessary, for example:
% This constructs the path for the above recording and the subfolder
% 'mousecam'
plab.locations.filename('server','AP001','2000-01-01','1200','mousecam')
% This constructs the path for the above animal and day in the folder
% 'ephys' (not in a recording time folder):
plab.locations.filename('server','AP001','2000-01-01',[],'ephys')
% And this example constructs a path 'histology' in the animal folder (not
% in a day folder): 
plab.locations.filename('server','AP001',[],[],'histology')

% Note that paths with plab.locations.filename are constructed whether or
% not the folder/file exists (e.g. above example paths do not exist).

% [EXERCISE] 
% Use 'plab.find_recording' to find the recording time for AP010 on
% 2023-08-10 with workflow 'lcr_passive', then use
% 'plab.locations.filename' to construct the path to the 'widefield'
% subfolder within the folder for that recording. Write a line to check
% whether that path exists on the server.

%% Loading data, and data types

% There isn't standardized lab code for loading data, since this is often
% customized for each person depending on their needs. This block demos my
% code, which can be used as a template, or as-is if it works for you (note
% that it's subject to regular changes).

% My loading code is in my repository (petersaj/AP_scripts_peterslab) and
% is called with 'ap.load_recording' after defining animal/day/time. This
% loads an example dataset with passive stimuli, widefield, and ephys. Run
% this data, then each data type will be explained below:
animal = 'AP005';
rec_day = '2023-06-21';
rec_time = '1859';
verbose = true; % this turns on/off progress display in command line
ap.load_recording;

% --- Timelite

% "Timelite" is our GUI for recording analog signals with a DAQ. These
% include things like wheel movement, stimulus screen photodiode, camera
% frame captures, etc. It is saved as:
% server\Data\animal\day\time\timelite.mat
% as the structure 'timelite':
timelite

% The structure 'timelite' contains these structures: 
% - daq_info: information about the DAQ recording
% - data: the recorded data as N timepoints x M signals
% - timestamps: the timestamps for each data point as N timepoints x 1

% The names of the recorded signals are under
% 'timelite.daq_info(channel number).channel_name'. For example, the
% photodiode is recorded under channel 5, and the name can be found by:
photodiode_channel = 5;
timelite.daq_info(photodiode_channel).channel_name

% The corresponding data is in the 5th column of 'data'. For example, this
% plots photodiode data by time: 
figure;plot(timelite.data(:,photodiode_channel));
xlabel('Time (s)');
ylabel('Voltage');
title('Photodiode signal');

% (The photodiode is over a square on the screen that turns white/black to
% indicate changes to the stimulus).

% [EXERCISE] 
% 1) Write code to identify the channel number for 'widefield_camera'
% (widefield camera exposures) and plot the data for that channel.
% 2) Find the timestamps when the widefield camera starts each exposure
% (times when the signal flips from low to high). 
% 3) Check that the timestamps in 2 match the variable
% 'widefield_expose_times' (this is created in ap.load_recording)

% --- Bonsai

% We use the program Bonsai to run our stimuli and tasks
% (https://bonsai-rx.org/). Bonsai is an actively supported, highly
% customizable, visual interfaced framework for designing experiments. 

% Bonsai files are called 'workflows'. ap.load_recording loads the name of
% the currently loaded workflow as 'bonsai_workflow'. This one is
% 'lcr_passive', which is passive stimuli presented on the left, center,
% and right screens.
bonsai_workflow

% Bonsai can change how it saves data, and loading scripts can specify how
% that data is loaded and parsed. In this demo, ap.load_recording loads
% data from Bonsai as 'trial_events':
trial_events
% which is a structure containing:
% - parameters: parameter values corresponding to the whole workflow, e.g.
% the duration that each stimulus was displayed was:
trial_events.parameters.StimDuration
% - values: an N trials x 1 array of values of saved events, e.g. the
% order of X-positions for the 3 stimuli presented on trial 4 was:
trial_events.values(4).TrialStimX
% - timestamps: an N trials x 1 array of timestamps of saved events from
% 'values', e.g. the time each of 3 stimuli was turned on/off on trial 4
% was:
trial_events.timestamps(4).StimOn

% Note that timing information in Bonsai is only approximate (when the
% software gave the command), not actual (when the stimulus was physically
% drawn on the screen). All time information should be used from Timelite,
% and only qualitative information should be used from Bonsai (e.g. which
% stimulus was presented on a trial).

% [EXERCISE] 
% The stimulus onset times are loaded by ap.load_recording as
% 'stimOn_times'. Write code to pull out a subset of these timestamps
% corresponding to a stimulus X position (TrialStimX) of 90 (meaning it was
% on the right-hand screen).

% --- Mousecam

% We record video of the front of the mouse during all experiments. The
% filename is loaded in ap.load_recording as:
mousecam_fn
% and the timestamps of each frame (in Timelite clock) is:
mousecam_times

% Mousecam images can be read into Matlab with a VideoReader object. For
% example:
mousecam_vr = VideoReader(mousecam_fn); % create VideoReader object

load_frame = 1; % define frame to read
mousecam_im = read(mousecam_vr,load_frame); % read frame

figure; % create a figure
imagesc(mousecam_im); % draw the image (sc = scaled colors)
axis image % make the x/y axes have equal aspect ratios
colormap('gray'); % set the colormap to gray

% You can also load in multiple frames at the same time by defining the
% start/end frame. For example:
load_frames = [1,10]; % define frame interval to read
mousecam_im = read(mousecam_vr,load_frames); % read frame
% (the VideoReader by default loads in multiple images in the 4th
% dimension, allowing for colors in the 3rd dimension. Our images are
% grayscale, so we can use 'squeeze' to remove the singleton 3rd dimension.
% Look at the size of the natively loaded data: 
size(mousecam_im)
% and then the size of the 'squeezed' data:)
mousecam_im = squeeze(mousecam_im); 
size(mousecam_im)
% This is an example of how to view a 3D matrix using my function
% ap.imscroll, which plots each 'page' (dim 1+2) which can be scrolled
% through in dim 3:
ap.imscroll(mousecam_im);

% [EXERCISE] 
% Create an average mousecam movie for -20:+20 frames around stimulus X =
% 90 presentations, and separately for stimulus X = -90. Is there a
% difference in behavior when these two stimuli are presented?

% --- Widefield and electrophysiology

% If available in a recording, ap.load_recording will load and prepare
% widefield and electrophysiology data. For demos on working with that
% data, see: 
% - plab_widefield_demo
% - plab_ephys_demo

% I have a function to scroll through experiment data for exploratory
% purposes, which displays the mousecam (left), widefield (center), and
% multiunit ephys (right), which can be scrolled through time:
ap.expscroll

% [EXERCISE] 
% There is a drop of sucrose available to the mouse at the beginning of
% this recording (left over from the task done previously). Which part of
% the ephys probe has increased activity when the mouse consumes this?

% --- Loading partial datasets

% Recordings can be part-loaded by ap.load_recording if only some data is
% necessary by creating a 'load_parts' structure, which can toggle loading
% with fields 'widefield','ephys','mousecam'. Timelite and behavior are
% always loaded. If 'load_parts' exists, any field not defined will default
% to not loading. For example, in this dataset, you can turn off loading of
% widefield and ephys data by doing:
clear all % (clears data loaded above)
animal = 'AP005';
rec_day = '2023-06-21';
rec_time = '1859';
load_parts.mousecam = true; % (this is the new line)
verbose = true;
ap.load_recording; % (this is much faster now and omits widefield and ephys)





















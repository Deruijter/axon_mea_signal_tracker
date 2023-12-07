%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code contains some snippets which should inform you how to use
%%% the AxtracktorClass
%%%
%%% Author: Markus de Ruijter; 2020; Punga Lab; Uppsala University




%% 1. Create AxtracktorClass instance (ALWAYS DO THIS (once)) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Add the location for all our libraries to the path
addpath(genpath('/mnt/nas/ClinicalNeurophysiology/software/matlab'));
% Set path to where the recording data is
path_data = '/mnt/argos/Silver/INV_ClinicalNeurophysiology/data/hdmea/';
% Create AxtracktorClass instance
axtracktor = AxtracktorClass(path_data);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NEW TRACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. Axtrackt all data that has not been axtrackted yet
% Note: this is what I (Markus) currently use to extract a new bunch of recordings

% Get path of all recordings
recording_id_list = axtracktor.GetRecordingIdListAll();
% Get "axdb" (table with all tracking data)
axdb = axtracktor.GetAxdb();

for i_ri = 1:size(recording_id_list,1)
	
	% Only axtrackt neurons that have not been manually checked/corrected yet
	if ~isempty(axdb(axdb.recording_id == recording_id_list.recording_id(i_ri) & axdb.track_quality ~= "unknown",:))
		continue;
	end
	
	% Create a new track (this will overwrite any previous data, even if it was manually corrected!!)
	axtracktor.NewTrack(recording_id_list.recording_id(i_ri));
	axtracktor.Axtrackt();
	
	disp('  ');
end
disp('done')





%% 3. Axtrackt specific recordings

% Set parameters of recordings you want to have tracked.
% An empty [] parameter indicates using ALL
% E.G. if you only set the date and the rest is empty, then axtrackt all recordings on that date
% If you set date, chip and experiment number => axtrackt all recordings for that date/chip/experiment
date = "2020_02_01";
chip_nr = 5195;
exp_nr = 2;
recording_nr = [1];	% [] = All 

% (optional) Check which recordings will be axtrackted:
recording_id_list = axtracktor.GetRecordingIdListAllWith(date, chip_nr, exp_nr, recording_nr)
% Perform axtracktion
axtracktor.AxtracktAllWith(date, chip_nr, exp_nr, recording_nr);	% This method overwrites existing tracks!!
disp('done')

% NOTE: Recordings that are in axdb, and do NOT have the track_quality "unknown" will NOT be axtrackted,
% or in other words, any recording that IS in axdb and has a track_quality other than "unknown" (e.g. "good"), will NOT be axtrackted






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXISTING TRACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4. Load an existing track

% Get the recording id for the neuron you want to load
date = "2020_02_01";
chip_nr = 5195;
exp_nr = 2;
recording_nr = 2;

recording_id_list = axtracktor.GetRecordingIdListAllWith(date, chip_nr, exp_nr, recording_nr)

axtracktor.LoadTrack(recording_id_list.recording_id);
disp('done')





%% 5. Re-axtrackt an existing track

% !!! CAREFUL, THIS WILL OVERWRITE ANY OLD DATA, INCLUDING MANUALLY CORRECTED DATA !!!

% Get the recording id for the neuron you want to re-axtrackt
date = "2020_02_01";
chip_nr = 5195;
exp_nr = 2;
recording_nr = 2;

recording_id_list = axtracktor.GetRecordingIdListAllWith(date, chip_nr, exp_nr, recording_nr)

axtracktor.NewTrack(recording_id_list.recording_id); % Using NewTrack() overwrites any old tracking data

% CHANGE SETTINGS HERE IF YOU WANT, SEE "8. Change various settings"

axtracktor.Axtrackt();
disp('done')





%% 6. View details of existing track

% Load existing track (See "4. Load an existing track")

% View peaks
peaks = axtracktor.GetPeaks();

% View connections/neurites
neurites = axtracktor.GetNeurites();

% View different recording formats
% "Raw" data (i.e. averaged voltage over multiple trials)
dt_voltage = axtracktor.GetDtVoltage();
% Amplitude data (converted from dt_voltage)
dt_amplitude = axtracktor.GetDtAmplitude();
% Tracking data (converted from dt_voltage, slightly different than dt_amplitude)
dt_tracking = axtracktor.GetDtTracking();

% View the axtracktor settings that were used during tracking
settings = axtracktor.GetSettings();

% View neuron statistics (same as what is stored in axdb)
neuron_stats = axtracktor.GetNeuronStats();
disp('done')





%% 7. Load axdb
% axdb (axtracktor database), contains all aggregated info of all tracked neurons

% No need to load any recording first

axdb = axtracktor.GetAxdb();






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 8. Change various settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See "5. Re-axtrackt an existing track" on how to incorporate these into a new axtracktion

%%
% Change the noise threshold using a modifier 
% 1 = Default threshold
% 0.5 = Half the noise threshold (= more peaks)
% 2 = Double the noise threshold (= less peaks)
% Can also be set to other values as long as they are above 0
axtracktor.SetSettingsNoiseThresholdModifier(1);

%%
% "Disable" advanced noise threshold
% This uses only one threshold (6 standard deviations).
% This was used to compare how well the advanced noise threshold worked VS a "standard" 6 standard deviation threshold
% Can theoretically be combined with SetSettingsNoiseThresholdModifier (don't think it makes much sense though)
% 1 = Use 6 standard deviation noise threshold (i.e. disable the advanced noise threshold)
% 0 = Use the advanced noise threshold (default)
axtracktor.SetSettingsUsePlainThreshTest(0);

%%
% Disable noise filter
% Used to assess the tracking with/without the noise filter
% 1 = use the noise filter (default)
% 0 = don't use the noise filter
axtracktor.SetSettingsNoiseFilterEnabled(1);

%%
% Set time span for amplitude calculation 
% The amplitude is the difference between the strongest and the weakest signal
% Here we can set what the time span is in which to look for the strongest/weakest signal
% 7 = default
axtracktor.SetSettingsDtAmplitudeTimeSpan(7);

%%
% Export graphics
% We can re-export the graphics, which will overwrite the old graphics files
% This takes a minute to finish
%"Time 1", [], 1,0,0,0
color_scheme = "Time 1"; % "Time 1" (default), "Speed", "Amplitude", "Time 2", "None" (don't draw lines)
alt_time_range = []; % [] = use all time frames; Specify a range to only output those time frames, e.g. [10:40]
make_eps = 1; % Boolean, output a picture of the final track (eps + png). Default = 1;
make_gif = 0; % Boolean, output an animated gif as well as a movie (can be useful for presentations). Default = 0;
plot_peaks = 0; % Boolean, output peaks as well as lines. Default = 0;
plot_special_peaks = 0; % Boolean, output "special" peaks (branches, end point, start point, initiation site). Default = 0;
axtracktor.ExportGraphics(color_scheme, alt_time_range, make_eps, make_gif, plot_peaks, plot_special_peaks);


%%
% Get a "span" of the movie file data
% This is used to view a certain "span" of the signal, e.g. the strongest signal over X time frames
% Must be combined with "ReloadDtAmplitude" or "ReloadDtTracking"
% DO NOT USE THIS WHEN TRACKING !!! ONLY FOR VISUALISATION
% 0 = don't use any time window (default)
% >0 = use X time frames to get the strongest signal over
axtracktor.SetSettingsDtTimeWindow(0);
% Example use:
recording_id_list = axtracktor.GetRecordingIdListAllWith("2020_02_01", 5195, 2, 2)
axtracktor.LoadTrack(recording_id_list.recording_id); % See "4. Load an existing track"
axtracktor.SetSettingsDtTimeWindow(30);
axtracktor.ReloadDtTracking(); % or ReloadDtAmplitude();
dt_tracking_tmsp10 = axtracktor.GetDtTracking();
imagesc(dt_tracking_tmsp10(:,:,70)); daspect([1,1,1]);	% you might want to log-scale it and flip it
% NOTE: This example is safe to use and does not override any data



























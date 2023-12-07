% Usage
% Make instance of the class:
% axtracktor = AxtracktorClass([path_to_data]);
%
% Get recording IDs
% recording_id_list = axtracktor.GetRecordingIdList();
%
% Create a new track for a recording:
% axtracktor.NewTrack(recording_id_list.recording_id(1))
% Or load an existing track:
% axtracktor.LoadTrack(recording_id_list.recording_id(1))
%
% Use public methods for processing, for example:
% axtracktor.Axtrackt(); % For doing complete automated axtracktion
% or
% peaks = axtracktor.GetPeaks(); % Get a table with all the peaks
% axtracktor.DelPeaks([1,2,3]); % Remove peaks with these IDs
%
% Author: Markus de Ruijter; 2020; Punga lab; Uppsala University

classdef AxtracktorClass < handle % See https://stackoverflow.com/questions/60759970/program-in-same-oop-style-as-app-designer

	%% ----- CONSTANTS -----
	properties (Constant)
		axtracktor_version = '2021.03.12';	% Only update for changes that affect the tracking. Format: yyyy.mm.dd
		neuron_types = ["motor", "cortical", "unknown"];
% 		track_qualities_name_full = ["unfixable_mult_neur", "unfixable_no_signal", "unfixable_noise", "fixable_hard", "fixable_easy", "complicated", "good", "unknown"];
% 		track_qualities_name_short = ["unf_mn", "unf_ns", "unf_n", "fix_h", "fix_e", "comp", "g", "unkn"];		
		track_qualities = table(["unfixable_mult_neur", "unfixable_no_signal", "unfixable_noise", "fixable_hard", "fixable_easy", "complicated", "good", "unknown"]' ...
			,["unf_mn", "unf_ns", "unf_n", "fix_h", "fix_e", "comp", "g", "unkn"]' ...
			,'variablenames',{'track_quality','track_quality_short'})
		
		% Noise thresholds. Standard deviations of the noise, as learned from the noise filtering test
		noise_thresholds_motor = struct('low',3,'medium',3.5,'high',8); 
		noise_thresholds_cortical = struct('low',3,'medium',3.5,'high',8);
	end

	%% ----- PROPERTIES (PRIVATE) -----
	properties (Access = private)
		axdb;
		path_data;
		recording_id;
		recording_id_list;
		remote_paths;
		
		% Neuron info
		neuron_type;
		
		% Settings
		settings = struct();
		
		% Data
		dt_voltage;
		dt_amplitude;
		dt_tracking;
		peaks;
		signal_init;
		neurites;
	end
	
		
	%% ----- GETTERS -----
	methods
		function ret = GetAxdb(obj)	
			if isempty(obj.axdb)
				obj.Output("GetAxdb()", 0);

				if isfile(strcat(obj.path_data, "axdb.csv"))
					obj.axdb = readtable(strcat(obj.path_data, "axdb.csv"), 'Delimiter','\t', 'TextType', 'string');

				else % Create nerdb if it doesn't exist
					variable_names_types = [["recording_id", "string"]; ...
											["axtracktor_version", "string"];...
											["track_quality", "string"]; ...
											["subtype", "string"]; ...
											["amplitude_initiation", "double"]; ...
											["amplitude_mean", "double"]; ...
											["amplitude_std", "double"]; ...
											["amplitude_median", "double"]; ...
											["amplitude_iqr", "double"]; ...
											["amplitude_q1", "double"]; ...
											["amplitude_q3", "double"]; ...
											["amplitude_d1", "double"]; ...
											["amplitude_d2", "double"]; ...
											["amplitude_d3", "double"]; ...
											["amplitude_d4", "double"]; ...
											["amplitude_d6", "double"]; ...
											["amplitude_d7", "double"]; ...
											["amplitude_d8", "double"]; ...
											["amplitude_d9", "double"]; ...
											["voltage_initiation", "double"]; ...
											["voltage_mean", "double"]; ...
											["voltage_std", "double"]; ...
											["voltage_median", "double"]; ...
											["voltage_iqr", "double"]; ...
											["voltage_q1", "double"]; ...
											["voltage_q3", "double"]; ...
											["voltage_d1", "double"]; ...
											["voltage_d2", "double"]; ...
											["voltage_d3", "double"]; ...
											["voltage_d4", "double"]; ...
											["voltage_d6", "double"]; ...
											["voltage_d7", "double"]; ...
											["voltage_d8", "double"]; ...
											["voltage_d9", "double"]; ...
											["speed_mean", "double"]; ...
											["speed_std", "double"]; ...
											["speed_median", "double"]; ...
											["speed_iqr", "double"]; ...
											["speed_q1", "double"]; ...
											["speed_q3", "double"]; ...
											["speed_d1", "double"]; ...
											["speed_d2", "double"]; ...
											["speed_d3", "double"]; ...
											["speed_d4", "double"]; ...
											["speed_d6", "double"]; ...
											["speed_d7", "double"]; ...
											["speed_d8", "double"]; ...
											["speed_d9", "double"]; ...
											["nr_endpoints", "double"]; ...
											["nr_branching_points", "double"]; ...
											["total_length", "double"]; ...
											["length_to_endpoints_mean", "double"]; ...
											["length_to_endpoints_std", "double"]; ...
											["length_to_endpoints_median", "double"]; ...
											["length_to_endpoints_iqr", "double"]; ...
											["length_to_endpoints_q1", "double"]; ...
											["length_to_endpoints_q3", "double"]; ...
											["inter_branch_length_mean", "double"]; ...
											["inter_branch_length_std", "double"]; ...
											["inter_branch_length_median", "double"]; ...
											["inter_branch_length_iqr", "double"]; ...
											["inter_branch_length_q1", "double"]; ...
											["inter_branch_length_q3", "double"]; ...
											["time_to_endpoints_mean", "double"]; ...
											["time_to_endpoints_std", "double"]; ...
											["time_to_endpoints_median", "double"]; ...
											["time_to_endpoints_iqr", "double"]; ...
											["time_to_endpoints_q1", "double"]; ...
											["time_to_endpoints_q3", "double"]]; 

					obj.axdb = table('Size',[0,size(variable_names_types,1)],... 
						'VariableNames', variable_names_types(:,1),...
						'VariableTypes', variable_names_types(:,2));
				end
			end
			ret = obj.axdb;
		end
		
		function ret = GetDtAmplitude(obj)
			if isempty(obj.dt_amplitude)
% 				disp('amp - a');
				obj.Output("GetDtAmplitude()", 0);
				
				tmsp = obj.GetSettings().dt_amplitude_time_span;
				tmwd = obj.GetSettings().dt_time_window; 
				dt_voltage = obj.GetDtVoltage(); % Create a local object for readability
% 				disp('amp - b');

				%tmsp = 1; % nr of time frames before and after the selected timeframe
				max_time = size(dt_voltage,3);
				dt_amplitude = zeros([size(dt_voltage,1), size(dt_voltage,2), min(size(dt_voltage,3)-tmsp,1)]);
% 				disp('amp - c');
				% Got to do this for each electrode (is there a faster way?)
				for elect_y = 1:size(dt_voltage,1)
					for elect_x = 1:size(dt_voltage,2)

						s1 = dt_voltage(elect_y, elect_x, :);

						% Get the min / max in a x timeframe timespan
						for i = 1:size(dt_voltage,3)
							time_window = max(i-tmwd,1):min(i+tmwd, max_time);
							[min_val, i_min_val] = min(s1(max(i-tmwd,1):min(i+tmwd, max_time)));
							i_min_time = time_window(i_min_val);

							% Slope before time i_min_time
							dt_amplitude(elect_y, elect_x, i) = max(  s1(  max(i_min_time-tmsp,1)  :  max(i_min_time-1, 1)  )  ) - min_val;
							% Slope after time i_min_time
							%amps(elect_y, elect_x, i) = (max(s1(min(i_min_time+1, max_time):min(i_min_time+tmsp, max_time))) - min_val);
							% Slope before + after i_min_time
% 							obj.dt_amplitude(elect_y, elect_x, i) = (max(s1(min(i_min_time+1, max_time):min(i_min_time+tmsp, max_time))) - min_val) + (max(s1(max(i_min_time-tmsp, 1):max(i_min_time-1, 1))) - min_val);
						end
					end
% 					disp(strcat("amp - d: ", num2str(elect_y)));
				end
				obj.dt_amplitude = dt_amplitude;
% 				disp('amp - z');
			end
			
			ret = obj.dt_amplitude;
		end
		
		function ret = GetDtTracking(obj)
			% Note, we calculate slopes here similar to in GetDtAmplitude() but slightly different. We're not interested in the amplitude but the most clear separation of signal vs noise.
			
			if isempty(obj.dt_tracking)
				obj.Output("GetDtTracking()", 99);
				
				tmsp = 1;
				tmwd = obj.GetSettings().dt_time_window;

				dt_voltage = obj.GetDtVoltage(); % Create a local object for readability. DtVoltage is basically the averaged raw data.

				%tmsp = 1; % nr of time frames before and after the selected timeframe
				max_time = size(dt_voltage,3);
				dt_tracking = zeros([size(dt_voltage,1), size(dt_voltage,2), min(size(dt_voltage,3)-tmsp,1)]);
				% Got to do this for each electrode (is there a faster way?)
				for elect_y = 1:size(dt_voltage,1)
% 						obj.Output(strcat(num2str(elect_y)), 99);
					for elect_x = 1:size(dt_voltage,2)

						s1 = dt_voltage(elect_y, elect_x, :);

						% Get the min / max in a x timeframe timespan
						for i_tf = 1:size(dt_voltage,3)

							time_window = min(max(i_tf-tmwd,1), max_time):max(min(i_tf+tmwd, max_time), 1); % If tmwd = 0, then this is just the current time frame
							[min_val, i_min_val] = min(s1(time_window)); % Get the lowest value in this time window (we want the lowest peak)
							i_min_time = time_window(i_min_val); % Time at which this minimal value was found

							% Slope before time i_min_time
							%obj.dt_amplitude(elect_y, elect_x, i) = (max(s1(max(i_min_time-tmsp, 1):max(i_min_time-1, 1))) - min_val);
							% Slope after time i_min_time
							%amps(elect_y, elect_x, i) = (max(s1(min(i_min_time+1, max_time):min(i_min_time+tmsp, max_time))) - min_val);
% 								if elect_y == 1 && elect_x == 1 & i_tf == 1 % output some info (for debugging)
% 									disp("DtTracking details:")
% 									reshape(s1, [size(dt_voltage,3),1])
% 									i_min_time
% 									max_time
% 									tmsp
% 									min_val
% 								end
							% Slope after + before i_min_time
							dt_tracking(elect_y, elect_x, i_tf) = (max(s1(min(i_min_time+1, max_time):min(i_min_time+tmsp, max_time))) - min_val) + (max(s1(max(i_min_time-tmsp, 1):max(i_min_time-1, 1))) - min_val);
						end
					end
				end
				obj.dt_tracking = dt_tracking;

				if obj.GetSettings().noise_filter.enabled
					if obj.GetNeuronType() == "cortical"
% 						obj.dt_tracking = obj.FilterBilateral(obj.dt_tracking);
% 						obj.dt_tracking = obj.FilterAdaptive(obj.dt_tracking);
						obj.dt_tracking = obj.FilterImgaussfilt3(obj.dt_tracking);
					else
						obj.dt_tracking = obj.FilterImgaussfilt3(obj.dt_tracking);
					end
				end
			end
			
			obj.dt_tracking(isnan(obj.dt_tracking)) = 0;
			
			ret = obj.dt_tracking;
		end
		
		function ret = GetDtVoltage(obj)
			if isempty(obj.dt_voltage)
				obj.Output("GetDtVoltage()", 0);
				
				% FRAMES
				if ~isfile(obj.GetRecordingId().path_struct) % We prefer to use the STRUCT instead of frames, this is the "old" way of reading the data
					% Get data file of averaged recordings with fixed positions
					mydata = load(obj.GetRecordingId().path_frames);

					% Shrink data    
					% For some reason the movie data files have different formats?
					if isfield(mydata, 'myMovie')
						mydata = mydata.myMovie.realtime;
					else
						mydata = mydata.myData;
					end

					dt_voltage = zeros([120, 220, size(mydata,2)]);
					for i = 1:size(mydata,2)
						d_cd = mydata{i}.CData;
						d_cd = imresize(d_cd, [120,220]);
						dt_voltage(:,:,i) = d_cd;
					end
					obj.dt_voltage = dt_voltage;
				end

				% STRUCT
				if isfile(obj.GetRecordingId().path_struct)
					data = load(obj.GetRecordingId().path_struct);
					data = data.myNeuron.temp;
					data(end:26400,:) = 0;	% For some reason we miss a data point some times? just set it to 0
					data = flip(rot90(reshape(data,[220,120,size(data,2)])),1);
					[a,peak_idx] = min(data,[],'all','linear');
					[a,b,peak_time] = ind2sub(size(data),peak_idx);

					min_time = max((peak_time - 20),1);
					max_time = min(peak_time + 100, size(data,3));
					data = data(:,:,min_time:max_time);

					obj.dt_voltage = data;
				end
			end
			
			ret = obj.dt_voltage;
		end
		
		function ret = GetNeurites(obj)
			if isempty(obj.neurites)
				obj.Output("GetNeurites()", 0);

				variable_names_types = [["x_from","double"];...
						["y_from","double"];...
						["x_to","double"];...
						["y_to","double"];...
						["time","double"];...
						["speed_el_tf","double"];...
						["length_el","double"];...
						["speed_ms","double"];...
						["amplitude","double"];...
						["voltage","double"];...
						["peak_id_from","double"];...
						["peak_id_to","double"];...
						["neurite_id","string"];...
						["segment_id","double"];...
						["branch_level","double"];...
						["cum_length","double"];...
						["is_end_point","double"];...
						["will_branch","double"];...
						["branch_length","double"];...
						["track_type","double"];...
						["previous_peak_time","double"]];
				obj.neurites = table('Size',[0,size(variable_names_types,1)],... 
					'VariableNames', variable_names_types(:,1),...
					'VariableTypes', variable_names_types(:,2));	
			end
				
			ret = obj.neurites;
		end
		
		function ret = GetNeuronStats(obj)
			% Note: not all details can be calculated (like track quality), these are stored in axdb
			obj.Output("GetNeuronStats()", 0);
			
			% Get data from axDB
			axdb = obj.GetAxdb();
			if sum(axdb.recording_id == obj.GetRecordingId().recording_id) == 0
				neuron_stats = obj.CreateEmptyTableCopy(obj.GetAxdb(), 1);
				neuron_stats.recording_id = obj.GetRecordingId().recording_id;
				obj.SetAxdbRecord(neuron_stats);	% Add a preliminary record to axdb
				
			elseif sum(axdb.recording_id == obj.GetRecordingId().recording_id) == 1
				neuron_stats = axdb(axdb.recording_id == obj.GetRecordingId().recording_id,:);
				
			elseif sum(axdb.recording_id == obj.GetRecordingId().recording_id) > 1
				error(strcat("ERROR: Multiple entries found in axdb for recording id: ", obj.GetRecordingId().recording_id, ...
					newline, "This shouldn't be possible. Remove one of the entries."));
			end
			
			neurites = obj.GetNeurites();
			peaks = obj.GetPeaks();
			signal_init = obj.GetSignalInit();
			settings = obj.GetSettings();
			
            neuron_subtypes = ["motor","cortical"];
            inter_branch_lengths = grpstats(neurites, {'neurite_id'},{'sum'},'datavars','length_el');
            
			% Calculate the statistics
            amplitude_quartiles = quantile(peaks.amplitude(peaks.time > (signal_init.time+8)), 3);
            amplitude_deciles = quantile(peaks.amplitude(peaks.time > (signal_init.time+8)), 9);
            voltage_quartiles = quantile(peaks.voltage(peaks.time > (signal_init.time+8)), 3);
            voltage_deciles = quantile(peaks.voltage(peaks.time > (signal_init.time+8)), 9);
            speed_quartiles = quantile(neurites.speed_ms, 3);
            speed_deciles = quantile(neurites.speed_ms, 9);
            length_to_endpoints_quartiles = quantile(neurites.cum_length(neurites.is_end_point==1)*17.5, 3);
            inter_branch_length_quartiles = quantile(inter_branch_lengths.sum_length_el*17.5, 3);
            time_to_endpoints_quartiles = quantile((neurites.time(neurites.is_end_point==1) - signal_init.time)*50, 3);
			
			if ismissing(neuron_stats.track_quality)	% If missing from axdb
				neuron_stats.track_quality = "unknown";
			end
			neuron_stats.axtracktor_version = convertCharsToStrings(settings.axtracktor_version);
			neuron_stats.subtype = obj.GetNeuronType(); 
            neuron_stats.amplitude_initiation = peaks.amplitude(peaks.x == signal_init.x & peaks.y == signal_init.y & peaks.time == signal_init.time);
            neuron_stats.amplitude_mean = mean(peaks.amplitude(peaks.time > (signal_init.time+8)));
            neuron_stats.amplitude_std = std(peaks.amplitude(peaks.time > (signal_init.time+8)));
            neuron_stats.amplitude_median = amplitude_quartiles(2);
            neuron_stats.amplitude_iqr = iqr(peaks.amplitude(peaks.time > (signal_init.time+8)));
            neuron_stats.amplitude_q1 = amplitude_quartiles(1);
            neuron_stats.amplitude_q3 = amplitude_quartiles(3);
            neuron_stats.amplitude_d1 = amplitude_deciles(1);
            neuron_stats.amplitude_d2 = amplitude_deciles(2);
            neuron_stats.amplitude_d3 = amplitude_deciles(3);
            neuron_stats.amplitude_d4 = amplitude_deciles(4);
            neuron_stats.amplitude_d6 = amplitude_deciles(6);
            neuron_stats.amplitude_d7 = amplitude_deciles(7);
            neuron_stats.amplitude_d8 = amplitude_deciles(8);
            neuron_stats.amplitude_d9 = amplitude_deciles(9);
            neuron_stats.voltage_initiation = peaks.voltage(peaks.x == signal_init.x & peaks.y == signal_init.y & peaks.time == signal_init.time);
            neuron_stats.voltage_mean = mean(peaks.voltage(peaks.time>(signal_init.time+8)));
            neuron_stats.voltage_std = std(peaks.voltage(peaks.time > (signal_init.time+8)));
            neuron_stats.voltage_median = voltage_quartiles(2);
            neuron_stats.voltage_iqr = iqr(peaks.voltage(peaks.time > (signal_init.time+8)));
            neuron_stats.voltage_q1 = voltage_quartiles(1);
            neuron_stats.voltage_q3 = voltage_quartiles(3);
            neuron_stats.voltage_d1 = voltage_deciles(1);
            neuron_stats.voltage_d2 = voltage_deciles(2);
            neuron_stats.voltage_d3 = voltage_deciles(3);
            neuron_stats.voltage_d4 = voltage_deciles(4);
            neuron_stats.voltage_d6 = voltage_deciles(6);
            neuron_stats.voltage_d7 = voltage_deciles(7);
            neuron_stats.voltage_d8 = voltage_deciles(8);
            neuron_stats.voltage_d9 = voltage_deciles(9);
            neuron_stats.speed_mean = mean(neurites.speed_ms);
            neuron_stats.speed_std = std(neurites.speed_ms);
            neuron_stats.speed_median = speed_quartiles(2);
            neuron_stats.speed_iqr = iqr(neurites.speed_ms);
            neuron_stats.speed_q1 = speed_quartiles(1);
            neuron_stats.speed_q3 = speed_quartiles(3);
            neuron_stats.speed_d1 = speed_deciles(1);
            neuron_stats.speed_d2 = speed_deciles(2);
            neuron_stats.speed_d3 = speed_deciles(3);
            neuron_stats.speed_d4 = speed_deciles(4);
            neuron_stats.speed_d6 = speed_deciles(6);
            neuron_stats.speed_d7 = speed_deciles(7);
            neuron_stats.speed_d8 = speed_deciles(8);
            neuron_stats.speed_d9 = speed_deciles(9);
            neuron_stats.nr_endpoints = sum(neurites.is_end_point);
            neuron_stats.nr_branching_points = sum(neurites.will_branch);
            neuron_stats.total_length = sum(neurites.length_el)*17.5;
            neuron_stats.length_to_endpoints_mean = mean(neurites.cum_length(neurites.is_end_point==1)*17.5);
            neuron_stats.length_to_endpoints_std = std(neurites.cum_length(neurites.is_end_point==1)*17.5);
            neuron_stats.length_to_endpoints_median = length_to_endpoints_quartiles(2);
            neuron_stats.length_to_endpoints_iqr = iqr(neurites.cum_length(neurites.is_end_point==1)*17.5);
            neuron_stats.length_to_endpoints_q1 = length_to_endpoints_quartiles(1);
            neuron_stats.length_to_endpoints_q3 = length_to_endpoints_quartiles(3);
            neuron_stats.inter_branch_length_mean = mean(inter_branch_lengths.sum_length_el*17.5);
            neuron_stats.inter_branch_length_std = std(inter_branch_lengths.sum_length_el*17.5);
            neuron_stats.inter_branch_length_median = inter_branch_length_quartiles(2);
            neuron_stats.inter_branch_length_iqr = iqr(inter_branch_lengths.sum_length_el*17.5);
            neuron_stats.inter_branch_length_q1 = inter_branch_length_quartiles(1);
            neuron_stats.inter_branch_length_q3 = inter_branch_length_quartiles(3);
            neuron_stats.time_to_endpoints_mean = mean((neurites.time(neurites.is_end_point==1) - signal_init.time)*50);
            neuron_stats.time_to_endpoints_std = std((neurites.time(neurites.is_end_point==1) - signal_init.time)*50);
            neuron_stats.time_to_endpoints_median = time_to_endpoints_quartiles(2);
            neuron_stats.time_to_endpoints_iqr = iqr((neurites.time(neurites.is_end_point==1) - signal_init.time)*50);
            neuron_stats.time_to_endpoints_q1 = time_to_endpoints_quartiles(1);
            neuron_stats.time_to_endpoints_q3 = time_to_endpoints_quartiles(3);
			
			ret = neuron_stats;
		end
		
		function ret = GetNeuronType(obj)
			if isempty(obj.neuron_type)
				obj.Output("GetNeuronType()", 0);

				% Get neuron type (is there a better way for this?)
				% Reading the 3rd line from the matlab script file, can't find any other
				% place this info is stored
				try
					matlab_script_files = dir(obj.GetRemotePaths().matlab);
					matlab_script_files = matlab_script_files(endsWith(lower({matlab_script_files.name}), '.m'),:);
					matlab_script_files = matlab_script_files(end); % Assume sorted by date, get the newest one
					fid=fopen(strcat(obj.remote_paths.matlab, matlab_script_files.name)); 
                    
					% set linenum to the desired line number that you want to import
					linenum = 3;
					% use '%s' if you want to read in the entire line or use '%f' if you want to read only the first numeric value
					linetext = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
					linetext = linetext{1}{1};
					linetext = regexp(linetext, '\w*', 'match');
					linetext = lower(strjoin(linetext));

					if contains(linetext, "motor")
						obj.neuron_type = obj.neuron_types(1);
					elseif contains(linetext, "cort")
						obj.neuron_type = obj.neuron_types(2);
					else
						obj.neuron_type = obj.neuron_types(3);
						obj.Output("Celltype could not be determined from script", 99);
					end
                catch
					obj.neuron_type = obj.neuron_types(3);
					obj.Output("Celltype could not be determined (no valid script found)", 99);
				end
			end
			ret = obj.neuron_type;
		end
		
		function ret = GetPeaks(obj)
			
			if isempty(obj.peaks) 
				obj.Output("GetPeaks()", 0);
				
				variable_names_types = [["id", "double"]; ...
										["x", "double"]; ...
										["y", "double"]; ...
										["time", "double"]; ...
										["voltage", "double"]; ...
										["amplitude", "double"]; ...
										["confidence", "double"]; ...
										["is_signal_init", "double"]];
				obj.peaks = table('Size',[0,size(variable_names_types,1)],... 
					'VariableNames', variable_names_types(:,1),...
					'VariableTypes', variable_names_types(:,2));
			end
			ret = obj.peaks;
		end
		
		function ret = GetRecordingId(obj)
			% Currently selected recording ID
			ret = obj.recording_id;
		end
		
		function ret = GetRecordingIdListAll(obj)
			% Iterate through Milos' file structure and get all recording id's
			% File structure is as follows:
			% [remote_server]/ClinicalNeurophysiology/milosh/[recording_date]/[chip_id]/spontscan/fixed/[experiment_number]/analyzes/movies/frames/*.mat
			% Author: Markus de Ruijter 2020
			%
			% Input
			%  pass_base_data: Location of the recordings (e.g. /mnt/nas/ClinicalNeurophysiology/milosh)

			path_chip_ids = '/';
			path_experiment_nr = '/spontscan/fixed/';
			path_neuron_nr = '/analyzes/struct/';

			variable_names = {'recording_id','recording_date','chip_id','experiment_nr','neuron_nr', 'file_prefix', 'path_tracking', 'path_tracking_grph', 'path_frames', 'path_movie', 'path_struct'};
			variable_types = {'string','string','double','double','double','string','string', 'string', 'string','string','string'};
			obj.recording_id_list = table('Size',[0,size(variable_names,2)],... % Table columns are explained at bottom of script
				'VariableNames', variable_names,...
				'VariableTypes', variable_types);


			recording_dates = dir(strcat(obj.path_data));
			recording_dates = recording_dates(~ismember({recording_dates.name},{'.','..'}));
			recording_dates = recording_dates([recording_dates.isdir]);
			% recording_dates

			for i_rd = 1:size(recording_dates,1)

				% Get chip_ids
				chip_ids = dir(strcat(recording_dates(i_rd).folder,'/', recording_dates(i_rd).name, path_chip_ids));
				chip_ids = chip_ids(~ismember({chip_ids.name},{'.','..'}));
				chip_ids = chip_ids([chip_ids.isdir]);
				if isempty(chip_ids)
					continue;
				end

				for i_ci = 1:size(chip_ids,1)

					% Get experiment numbers
					experiment_nrs = dir(strcat(chip_ids(i_ci).folder,'/', chip_ids(i_ci).name, path_experiment_nr));
					experiment_nrs = experiment_nrs(~ismember({experiment_nrs.name},{'.','..'}));
					experiment_nrs = experiment_nrs([experiment_nrs.isdir]);
					if isempty(experiment_nrs)
						continue;
					end

					for i_en = 1:size(experiment_nrs,1)
						% Get experiment numbers
						neuron_nrs = dir(strcat(experiment_nrs(i_en).folder,'/', experiment_nrs(i_en).name, path_neuron_nr));
						neuron_nrs = neuron_nrs(~ismember({neuron_nrs.name},{'.','..'}));
						neuron_nrs = neuron_nrs(~[neuron_nrs.isdir]);
						if isempty(neuron_nrs)
							continue;
						end

						for i_ni = 1:size(neuron_nrs,1)
														
							chip_id = join(regexp(chip_ids(i_ci).name, ['\d'], 'match'),'');
							chip_id = str2num(chip_id{1});
					
							exp_nr = join(regexp(experiment_nrs(i_en).name, ['\d'], 'match'),'');
							exp_nr = str2num(exp_nr{1});
					
							neuron_nr = join(regexp(neuron_nrs(i_ni).name, ['\d'], 'match'),'');
							neuron_nr = str2num(neuron_nr{1});
							
							recording_id = strcat('d',recording_dates(i_rd).name,'_c', num2str(chip_id),'_e', num2str(exp_nr,'%02.f'),'_r', num2str(neuron_nr,'%02.f'));
							
							file_prefix = regexp(neuron_nrs(i_ni).name,'[^\d\W]*','match'); % Get all the words (not numbers)
                            file_prefix = file_prefix(1); % First word we found is the file prefix
							path_fileshelf = fileshelf('dirname', recording_dates(i_rd).name, 'chipID', chip_id, 'spontscan', 'fixed',  'expNo', exp_nr, 'basepath', obj.path_data);
							path_tracking = strcat(path_fileshelf.axt_data, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
							path_tracking_grph = strcat(path_fileshelf.axt_grph, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
							path_frames = strcat(path_fileshelf.frames, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
							path_movie = strcat(path_fileshelf.movfiles, 'myMovie', num2str(neuron_nr,'%02.f'),'.avi');
							path_struct = strcat(path_fileshelf.struct, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
			
							obj.recording_id_list(end+1,:) = {recording_id, recording_dates(i_rd).name, chip_id, exp_nr,  neuron_nr, file_prefix, path_tracking, path_tracking_grph, path_frames, path_movie, path_struct};
						end
					end
				end
			end
			obj.recording_id_list = obj.recording_id_list(~contains(lower(obj.recording_id_list.file_prefix), 'muscle'),:);
			
			% Append axdb and track quality data
			obj.recording_id_list = outerjoin(obj.recording_id_list, obj.GetAxdb(), 'Type','left','Keys','recording_id','MergeKeys',1);
% 			obj.recording_id_list
% 			obj.recording_id_list = join(obj.recording_id_list, obj.track_qualities, 'LeftKeys','track_quality','RightKeys','track_quality');
			obj.recording_id_list = outerjoin(obj.recording_id_list, obj.track_qualities, 'Type','left','Keys','track_quality','MergeKeys',1);
			
			ret = obj.recording_id_list; % Skip muscle cells for now
		end
		
		function ret = GetRecordingIdListAllWith(obj, date, chip_nr, exp_nr, recording_nr)
			% Function specifically made to work with Milos' filepaths
			% We don't have chip/date/experiment info, only the path
			% Arguments may be empty in a hierarchical manner:
			% [],[],[],[]						=> Same as GetRecordingIdListAll()
			% '2020_11_11',[],[],[]				=> Get all recordings on this date
			% '2020_11_11',1234,[],[]			=> Get all recordings for this chip on this date
			% etc..

			path_chip_ids = '/';
			path_experiment_nr = '/spontscan/fixed/';
			path_neuron_nr = '/analyzes/struct/';

			variable_names = {'recording_id','recording_date','chip_id','experiment_nr','neuron_nr', 'file_prefix', 'path_tracking', 'path_tracking_grph', 'path_frames', 'path_movie', 'path_struct'};
			variable_types = {'string','string','double','double','double','string','string', 'string', 'string','string','string'};
			obj.recording_id_list = table('Size',[0,size(variable_names,2)],... % Table columns are explained at bottom of script
				'VariableNames', variable_names,...
				'VariableTypes', variable_types);


			if isempty(date)
				recording_dates = dir(strcat(obj.path_data));
				recording_dates = recording_dates(~ismember({recording_dates.name},{'.','..'}));
				recording_dates = recording_dates([recording_dates.isdir]);
			else
				recording_dates(1).name = date;
			end
			
			% recording_dates
			for i_rd = 1:size(recording_dates,1)
				cur_folder = strcat(obj.path_data, recording_dates.name(i_rd), path_chip_ids);

				if isempty(chip_nr)
					% Get chip_ids
					chip_ids = dir(cur_folder);
					chip_ids = chip_ids(~ismember({chip_ids.name},{'.','..'}));
					chip_ids = chip_ids([chip_ids.isdir]);
					if isempty(chip_ids)
						continue;
					end
				else
					chip_ids(1).name = strcat('chip',num2str(chip_nr));
				end

				for i_ci = 1:size(chip_ids,1)
					cur_folder = strcat(cur_folder, chip_ids(i_ci).name, path_experiment_nr);

					if isempty(exp_nr)
						% Get experiment numbers
						experiment_nrs = dir(cur_folder);
						experiment_nrs = experiment_nrs(~ismember({experiment_nrs.name},{'.','..'}));
						experiment_nrs = experiment_nrs([experiment_nrs.isdir]);
						if isempty(experiment_nrs)
							continue;
						end
					else
						experiment_nrs(1).name = strcat('experiment_',num2str(exp_nr,'%02.f'));
					end

					for i_en = 1:size(experiment_nrs,1)
						cur_folder = strcat(cur_folder, '/', experiment_nrs(i_en).name, path_neuron_nr);
						
						if isempty(recording_nr)
							% Get experiment numbers
							recording_nrs = dir(cur_folder);
							recording_nrs = recording_nrs(~ismember({recording_nrs.name},{'.','..'}));
							recording_nrs = recording_nrs(~[recording_nrs.isdir]);
							if isempty(recording_nrs)
								continue;
							end
							file_prefix = regexp(recording_nrs(1).name,'[^\d\W]*','match'); % Get all the words (not numbers)
                            file_prefix = file_prefix(1); % First word we found is the file prefix
						else
							% We need the file prefix...
							tmp = dir(cur_folder);
							tmp = tmp(~ismember({tmp.name},{'.','..'}));
							tmp = tmp(~[tmp.isdir]);
							if isempty(tmp)
								continue;
							end
							file_prefix = regexp(tmp(1).name,'[^\d\W]*','match'); % Get all the words (not numbers)
                            file_prefix = file_prefix{1}; % First word we found is the file prefix
							
							recording_nrs(1).name = strcat(file_prefix, num2str(recording_nr,'%02.f')); 							
						end

						for i_ni = 1:size(recording_nrs,1)
														
							chip_id = join(regexp(chip_ids(i_ci).name, ['\d'], 'match'),'');
							chip_id = str2num(chip_id{1});
					
							exp_nr = join(regexp(experiment_nrs(i_en).name, ['\d'], 'match'),'');
							exp_nr = str2num(exp_nr{1});
					
							neuron_nr = join(regexp(recording_nrs(i_ni).name, ['\d'], 'match'),'');
							neuron_nr = str2num(neuron_nr{1});
							
							recording_id = strcat('d',recording_dates(i_rd).name,'_c', num2str(chip_id),'_e', num2str(exp_nr,'%02.f'),'_r', num2str(neuron_nr,'%02.f'));
							
							path_fileshelf = fileshelf('dirname', convertStringsToChars(recording_dates(i_rd).name), 'chipID', chip_id, 'spontscan', 'fixed',  'expNo', exp_nr, 'basepath', obj.path_data);
							path_tracking = strcat(path_fileshelf.axt_data, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
							path_tracking_grph = strcat(path_fileshelf.axt_grph, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
							path_frames = strcat(path_fileshelf.frames, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
							path_movie = strcat(path_fileshelf.movfiles, 'myMovie', num2str(neuron_nr,'%02.f'),'.avi');
							path_struct = strcat(path_fileshelf.struct, file_prefix, num2str(neuron_nr,'%02.f'),'.mat');
			
							obj.recording_id_list(end+1,:) = {recording_id, recording_dates(i_rd).name, chip_id, exp_nr,  neuron_nr, file_prefix, path_tracking, path_tracking_grph, path_frames, path_movie, path_struct};
						end
					end
				end
			end
			obj.recording_id_list = obj.recording_id_list(~contains(lower(obj.recording_id_list.file_prefix), 'muscle'),:); % Skip muscle cells for now
			
			% Append axdb and track quality data
			obj.recording_id_list = outerjoin(obj.recording_id_list, obj.GetAxdb(), 'Type','left','Keys','recording_id','MergeKeys',1);
			obj.recording_id_list = outerjoin(obj.recording_id_list, obj.track_qualities, 'Type','left','Keys','track_quality','MergeKeys',1);
			
			ret = obj.recording_id_list; 
		end
				
		function ret = GetRemotePaths(obj)
			if isempty(obj.remote_paths)
				obj.Output("GetRemotePaths()", 0);

				recording_date = convertStringsToChars(obj.GetRecordingId().recording_date);
				chip_id = obj.GetRecordingId().chip_id;
				experiment_nr = obj.GetRecordingId().experiment_nr;
				obj.remote_paths = fileshelf('basepath', obj.path_data, 'dirname',recording_date,'chipID',chip_id,'fixed','expNo',experiment_nr); % External function, see fileshelf.m
			end
			ret = obj.remote_paths;
		end
	
		function ret = GetSettings(obj)
			
			if isempty(obj.settings) || isempty(fieldnames(obj.settings))  
				obj.Output("GetSettings()", 0);
			
				obj.settings.max_speed = 6; % electrodes per time unit. NOTE: Need to actually check this.
				obj.settings.peak_conn = 8; % conn parameter for imregionalmax, 4 = only perpendicular; 8 = also diagonally
				obj.settings.neuron_type = obj.neuron_type;
				obj.settings.axtracktor_version = obj.axtracktor_version;
				
				obj.settings.dt_amplitude_time_span = 7; % The time span over which to determine the amplitude
				obj.settings.dt_time_window = 0;	% Used to create a "span" of the video instead of measurements of individual time frames (0 = use no span, only individual frames)

				% Peak detection settings
				obj.settings.noise_threshold_modifier = 1;
				obj.settings.noise_threshold.motor.std.high = obj.noise_thresholds_motor.high;
				obj.settings.noise_threshold.motor.std.medium = obj.noise_thresholds_motor.medium;
				obj.settings.noise_threshold.motor.std.low = obj.noise_thresholds_motor.low;
				obj.settings.noise_threshold.cortical.std.high = obj.noise_thresholds_cortical.high;
				obj.settings.noise_threshold.cortical.std.medium = obj.noise_thresholds_cortical.medium;
				obj.settings.noise_threshold.cortical.std.low = obj.noise_thresholds_cortical.low;
				obj.settings.noise_threshold.unknown.std.high = obj.noise_thresholds_motor.high;
				obj.settings.noise_threshold.unknown.std.medium = obj.noise_thresholds_motor.medium;
				obj.settings.noise_threshold.unknown.std.low = obj.noise_thresholds_motor.low;
				obj.settings.noise_threshold.used.std.high = 0;
				obj.settings.noise_threshold.used.std.medium = 0;
				obj.settings.noise_threshold.used.std.low = 0;
				obj.settings.noise_threshold.used.absolute.high = 0;
				obj.settings.noise_threshold.used.absolute.medium = 0;
				obj.settings.noise_threshold.used.absolute.low = 0;
				obj.settings.noise_threshold.motor.absolute.high = 0;
				obj.settings.noise_threshold.motor.absolute.medium = 0;
				obj.settings.noise_threshold.motor.absolute.low = 0;
				obj.settings.noise_threshold.cortical.absolute.high = 0;
				obj.settings.noise_threshold.cortical.absolute.medium = 0;
				obj.settings.noise_threshold.cortical.absolute.low = 0;
% 				obj.settings.noise_threshold.unknown.absolute.high = 0;
% 				obj.settings.noise_threshold.unknown.absolute.medium = 0;
% 				obj.settings.noise_threshold.unknown.absolute.low = 0;
				obj.settings.neighbor_dist.high = 0;   
				obj.settings.neighbor_dist.medium = 3;
				obj.settings.neighbor_dist.low = 1;
				obj.settings.neighbor_lookback_time = 2;
				obj.settings.noise_threshold.use_plain_thresh_test = 0;	% For testing how much the other thresholds improve predictions
				obj.settings.noise_threshold.plain_thresh_test.std.high = 6;
				obj.settings.noise_threshold.plain_thresh_test.std.medium = inf;
				obj.settings.noise_threshold.plain_thresh_test.std.low = inf;

				% Noise filters
				obj.settings.noise_filter.enabled = 1;
				obj.settings.noise_filter.adaptive.size_1 = 11;
				obj.settings.noise_filter.adaptive.size_2 = 11;
				obj.settings.noise_filter.imgaussfilt3.sigma = 0.7;

				% Peak connection
				obj.settings.connect_look_back_time = 5; % time frames
				obj.settings.connect_direct_min_peak_dist = 3; % Direct connection distance (i.e. don't use dynamic programming)

				% Dynamic programming
				obj.settings.use_dynamic_programming = 1;
				obj.settings.dp.radius = 13;
				obj.settings.dp.time_range = 18;
				obj.settings.dp.rank_weights.mean = 1;
				obj.settings.dp.rank_weights.q_1 = 1;
				obj.settings.dp.rank_weights.std = 2;
				obj.settings.dp.rank_weights.speed_diff = 0.5;
				obj.settings.dp.rank_weights.distance = 4;
				obj.settings.dp.rank_weights.angle = 4;
			end
			
			ret = obj.settings;
		end
		
		function ret = GetSignalInit(obj)
			peaks = obj.GetPeaks();
			signal_init = peaks(peaks.is_signal_init == 1,:);
			
			if isempty(signal_init)
				obj.Output("GetSignalInit()", 0);
				
				dt_amplitude = obj.GetDtAmplitude();
				dt_voltage = obj.GetDtVoltage();

				signal_init = obj.CreateEmptyTableCopy(peaks,1);
% 				dt_tracking = obj.GetDtTracking();
				signal_start = 1;
				highest_signal_overall = 0;
				for i = 1:size(dt_amplitude,3)
					highest_signal = max(dt_amplitude(:,:,i),[],  'all');
					if highest_signal > highest_signal_overall
					   highest_signal_overall = highest_signal;
					   signal_start = i;
					end
				end
				dt_ss = dt_amplitude(:,:,signal_start); %myData{signal_start}.CData;
				[row, col] = find(dt_ss == max(dt_ss(:)));  % assume we have only one peak for now

				signal_init.x = round(median(col));
				signal_init.y = round(median(row));
				signal_init.time = signal_start;
				signal_init.amplitude = dt_amplitude(signal_init.y,signal_init.x,signal_start);
				signal_init.voltage = dt_voltage(signal_init.y,signal_init.x,signal_start);
				signal_init.confidence = 0; % TODO
				signal_init.is_signal_init = 1;
			end
			
			ret = signal_init;
		end
	end
	
	%% ----- SETTERS -----
	methods
		function SetAxdbRecord(obj, neuron_stats)
			obj.Output("SetAxdbRecord()", 0);
			if sum(obj.axdb.recording_id == obj.GetRecordingId().recording_id) == 0
				obj.axdb(end+1,:) = neuron_stats;
				
			elseif sum(obj.axdb.recording_id == obj.GetRecordingId().recording_id) == 1
				obj.axdb(obj.axdb.recording_id == obj.GetRecordingId().recording_id,:) = neuron_stats;
				
			else sum(obj.axdb.recording_id == obj.GetRecordingId().recording_id) > 1
				error(strcat("ERROR: Multiple entries found in axdb for recording id: ", obj.GetRecordingId().recording_id, ...
					newline, "This shouldn't be possible. Remove one of the entries."));
			end
		end
		
		function SetRecordingId(obj, recording_id)
			% NOTE: All variables are reset when using this method (including settings)
			obj.Output("SetRecordingId()", 0);
			
			recording_id_list = obj.recording_id_list;
			obj.recording_id = recording_id_list(recording_id_list.recording_id == recording_id,:);
			obj.settings = [];
			obj.peaks = [];
			obj.neurites = [];
			obj.dt_voltage = [];
			obj.dt_amplitude = [];
			obj.dt_tracking = [];
		end
		
		function SetSettingsDtTimeWindow(obj, time_window)
			% Set the time window for GetDtAmplitude and GetDtTracking
			% This should not be used when performing the tracking!
			% Only for viewing the data !!
			disp("Don't use SetSettingsDtTimeWindow() for tracking!")
			
			settings = obj.GetSettings();
			settings.dt_time_window = time_window;
			obj.settings = settings;
		end
	
		function SetSettingsNoiseFilterEnabled(obj, on)
			obj.Output("SetSettingsNoiseFilterOnOff()", 0);
			
			settings = obj.GetSettings();
			settings.noise_filter.enabled = on;
			obj.settings = settings;
		end
		
		function SetSettingsDisableComplexRouting(obj)
			% Reset to default settings
			settings = obj.GetSettings();
			settings.use_dynamic_programming = 0;
			obj.settings = settings;
		end
	
		function SetSettingsDtAmplitudeTimeSpan(obj, time_span)
			obj.Output("SetSettingsNoiseFilterOnOff()", 0);
			
			settings = obj.GetSettings();
			settings.dt_amplitude_time_span = time_span;
			obj.settings = settings;
		end
		
		function SetSettingsNoiseThresholdModifier(obj, modifier)
			obj.Output("SetNoiseThresholdModifier()", 0);
			
			% Increase or decrease the default noise threshold
			% with a certain percentage
			settings = obj.GetSettings();
			settings.noise_threshold_modifier = modifier;
			obj.settings = settings;			
		end
		
		function SetSettingsUsePlainThreshTest(obj, on)
			obj.Output("SetSettingsUseNoThreshTest()", 0);
			
			settings = obj.GetSettings();
			settings.noise_threshold.use_plain_thresh_test = on;
			obj.settings = settings;
		end
		
		function SetTrackQuality(obj, quality)
			obj.Output("SetTrackQuality()", 0);
			
			if ismember(quality, obj.track_qualities.track_quality)
				neuron_stats = obj.GetNeuronStats();
				neuron_stats.track_quality = quality;
				obj.SetAxdbRecord(neuron_stats);
			else
				disp("Error - Track quality must be one of following:")
				disp(strjoin(obj.track_qualities.name_full, " "));
			end
		end
		
		function SetRankedRouteWeights(obj, param, weight)
			obj.Output("SetRankedRouteWeight()", 0);
			
			
			settings = obj.GetSettings();
			settings.dp.rank_weights.(param) = weight;
			obj.settings = settings;
		end
		
		function SetSettingsDefault(obj)
			% Reset to default settings
			
			obj.settings = [];
			settings = obj.GetSettings();
		end
	end
	

	%% ----- METHODS (PUBLIC) -----
	methods
		% Constructor
		function obj = AxtracktorClass(path_data)
			% path_data = String, path to the recording files
			obj.path_data = path_data;
		end
		
		% Main function, perform axtracktion
		function ret = Axtrackt(obj)
			obj.Output("Axtrackt()", 0);
			obj.Output(strcat("Version ", obj.axtracktor_version), 0);
			
			if ~obj.IsValidRecordingId(); return; end
			if ~obj.HasValidCelltype(); return; end     % Don't process unknown cell types
			%if ~obj.IsToBeAxtrackted(); return; end
			
			% Perform peak finding & tracking
			obj.Output("Axtrackt() a", 99);
			obj.ClearTrackingData();
			obj.Output("Axtrackt() b", 99);
			obj.FindPeaks();
			obj.Output("Axtrackt() c", 99);
			obj.RemovePrePostPeaks();
			obj.Output("Axtrackt() d", 99);
			obj.ConnectPeaks();
			obj.Output("Axtrackt() e", 99);
			obj.AppendNeuriteLengthsAndIdsAndStuff();
			obj.Output("Axtrackt() f", 99);
			obj.RemoveStationaryPeaks();
			obj.Output("Axtrackt() g", 99);
			
			% Export / save
			obj.ExportGraphics("Time 1", [], 1,0,0,0); % <- Maybe split up into different methods
			obj.ExportData();
			obj.Output("Axtrackt() h", 99);
			
			% Add or update the statistics of this neuron in axon database
			obj.SetAxdbRecord(obj.CalculateNeuronStats());
			obj.SaveAxdb();
			
			ret = 1;
		end
		
		% Main function wrapper, select the date/chip/experiment/recording nr you want recorded
		% Leave one of these empty to axtrackt all of them (e.g. empty recording_nr => axtrackt all
		% recordings in that date/chip/experiment)
		function ret = AxtracktAllWith(obj, date, chip_nr, exp_nr, recording_nr)
			% Note: This function is made specifically to be used within Milos' "base script"
			
			recording_id_list = obj.GetRecordingIdListAllWith(date, chip_nr, exp_nr, recording_nr);
			
			axdb = obj.GetAxdb();
			
			for i_ri = 1:size(recording_id_list,1)
% 				axtracktor_tic = tic;

				disp(['Axtrackting ' num2str(i_ri) ' of ' num2str(size(recording_id_list,1)) ': ' convertStringsToChars(recording_id_list.recording_id(i_ri))]);

				% Only axtrackt neurons that have not been manually checked/corrected yet
				if ~isempty(axdb(axdb.recording_id == recording_id_list.recording_id(i_ri) & axdb.track_quality ~= "unknown",:))
					continue;
				end

				% Create a new track (this will overwrite any previous data, even if it was manually corrected!!)
				obj.NewTrack(recording_id_list.recording_id(i_ri));
				obj.Axtrackt();

% 				toc(axtracktor_tic);
			end
			
			ret = 1;
		end
		
		function nul = AddPeaks(obj, new_peaks)
			% new_peaks should be a table with the same columns as obj.GetPeaks()
			% use CreateEmptyTableCopy() to create an empty copy of a table
			if ~isempty(new_peaks)
				peaks = obj.GetPeaks();
				new_peaks = new_peaks(~ismember(new_peaks{:,{'x','y','time'}}, peaks{:,{'x','y','time'}},'rows'),:);
				max_id = max([obj.GetPeaks().id ; 0]);
				new_ids = (max_id+1) : (max_id+1)+(size(new_peaks,1)-1);
				new_peaks.id = new_ids';
				obj.peaks = [peaks ; new_peaks];
			end
		end
		
		function nul = ClearVariables(obj)
			obj.dt_voltage = [];
			obj.dt_amplitude = [];
			obj.dt_tracking = [];
			obj.neuron_type = [];
			obj.remote_paths = [];
			obj.settings = [];
			obj.peaks = [];
			obj.neurites = [];
			obj.signal_init = [];
		end
		
		function nul = DelPeaks(obj, peak_id_list)
			obj.Output("DelPeaks()", 0);
			
			if isempty(peak_id_list)
				return;
			end
			
			% Not allowed to remove the initiation site
			peaks = obj.GetPeaks();
			signal_init_peak_id = peaks.id(peaks.is_signal_init == 1);
			peak_id_list = peak_id_list(peak_id_list ~= signal_init_peak_id);
			
			obj.peaks(ismember(obj.peaks.id, peak_id_list),:) = []; 
			
			% Remove associated connections
			obj.DelNeuritesFrom(peak_id_list);
			obj.DelNeuritesTo(peak_id_list);
		end
		
		function ret = DelPeaksAndAfter(obj, peak_id_list)
			% Delete all peaks in the list, along with any later peak that is connected to this peak.
			% ret = list of deleted peak ids
			obj.Output("DelPeaksAndAfter()", 0);
			
			if isempty(peak_id_list)
				return;
			end
			
			% Not allowed to remove the initiation site
			neurites = obj.GetNeurites();
			peaks = obj.GetPeaks();
			
			branch_id = neurites(ismember(neurites.peak_id_to, peak_id_list), {'neurite_id','segment_id'});
			peak_ids_after = nan([99999,1]);	% Preallocate space to prevent memory leaks
			pia_idx = 1;
			for i_pi = 1:size(branch_id,1)
				peak_ids_after_tmp = neurites{(neurites.neurite_id == branch_id.neurite_id(i_pi) & neurites.segment_id >= branch_id.segment_id(i_pi)) | (startsWith(neurites.neurite_id, branch_id.neurite_id(i_pi)) & neurites.neurite_id ~= branch_id.neurite_id(i_pi)),'peak_id_to'}
				peak_ids_after(pia_idx:(pia_idx+size(peak_ids_after_tmp,1)-1)) = peak_ids_after_tmp';
				pia_idx = pia_idx+size(peak_ids_after_tmp)
			end
			peak_id_list = unique(peak_ids_after(~isnan(peak_ids_after)))

			signal_init_peak_id = peaks.id(peaks.is_signal_init == 1);
			peak_id_list = peak_id_list(peak_id_list ~= signal_init_peak_id);
			
			ret = obj.peaks{ismember(obj.peaks.id, peak_id_list),'id'};
			
			obj.peaks(ismember(obj.peaks.id, peak_id_list),:) = []; 
			
			% Remove associated connections
			obj.DelNeuritesFrom(peak_id_list);
			obj.DelNeuritesTo(peak_id_list);
		end
		
		function nul = AddNeurites(obj, new_neurites)
			% TODO: Set IDs
			obj.neurites = [obj.GetNeurites() ; new_neurites];
		end
		
		function nul = DelNeuritesTo(obj, peak_id_list)
			obj.Output("DelNeuritesTo()", 0);
			
			if isempty(peak_id_list) || isempty(obj.neurites)
				return;
			end
			
			obj.neurites(ismember(obj.neurites.peak_id_to, peak_id_list),:) = [];
			obj.AppendNeuriteLengthsAndIdsAndStuff();
		end
		
		function nul = DelNeuritesFrom(obj, peak_id_list)
			obj.Output("DelNeuritesFrom()", 0);
			
			if isempty(peak_id_list) || isempty(obj.neurites)
				return;
			end
			
			obj.neurites(ismember(obj.neurites.peak_id_from, peak_id_list),:) = [];
			obj.AppendNeuriteLengthsAndIdsAndStuff();
		end
		
		function nul = SaveAxdb(obj)
			obj.Output("SaveAxdb()", 0);
			
			axdb = obj.GetAxdb();
            writetable(axdb, strcat(obj.path_data,'axdb.csv'), 'Delimiter','\t');
		end
		
		function nul = ExportGraphics(obj, color_scheme, alt_time_range, make_eps, make_gif, plot_peaks, plot_special_peaks)
            % Note: Should probably split this method up into smaller chunks at some point
			obj.Output("ExportGraphics()", 0);
			
			neurites = obj.GetNeurites();
			peaks = obj.GetPeaks();
			settings = obj.GetSettings();
			signal_init = obj.GetSignalInit();
			dt_filtered = obj.GetDtTracking();
			remote_paths = obj.GetRemotePaths();
			recording_id = obj.GetRecordingId();
			
			%neurites = table2array(neurites);
			neurites = sortrows(neurites, 5);
			null = [];
			% Create avi file
			%disp('Creating video');
			% Use subtightplot 
			make_it_tight = true;
			subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.05], 0.025, 0.05);
			%subplot = @(m,n,p) subtightplot (m, n, p, [0 0], 0, 0);
			if ~make_it_tight,  clear subplot;  end

			electrodes_x = size(dt_filtered,2);
			electrodes_y = size(dt_filtered,1);
			time_points = size(dt_filtered,3);

			file_base_name = strcat(obj.GetRecordingId().path_tracking_grph);
			
			videoName = strcat(file_base_name, '.avi');
			gifName = strcat(file_base_name, '.gif');
			svgName = strcat(file_base_name, '.eps');
			pngName = strcat(file_base_name, '.png');
% 			loops = 10;
% 			F(loops) = struct('cdata',[],'colormap',[]);
			timediff = 1;
			if ~isempty(neurites)
				timediff = max(neurites.time - min(neurites.time)+1);
			end
			colcode.realtime = jet(timediff); %spring;
			background_cls = colcode.realtime;


			line_cls = colcode.realtime(neurites.time - min(neurites.time)+1,:);
% 			if color_scheme == 2
% 				line_cls = repmat([0,1,1], size(neurites,1), 1);
% 				line_cls(neurites.track_type == 2, 1) = 1.0;   % Dynamic programming connections
% 				line_cls(neurites.track_type == 2, 2) = 1.0;
% 				line_cls(neurites.track_type == 2, 3) = 1.0;
% 				line_cls(neurites.track_type == 3, 1) = 0.5; % Pseudo connections
% 				line_cls(neurites.track_type == 3, 2) = 0.5;
% 				line_cls(neurites.track_type == 3, 3) = 0.5;
% 				background_cls = line_cls;
% 
% 			elseif color_scheme == 3    % Better for presentations / projector
% 				colcode.realtime = spring(max(neurites.time - min(neurites.time)+1)); %spring;
% 				background_cls = colcode.realtime;
% 				line_cls = colcode.realtime(neurites.time - min(neurites.time)+1,:);
% 
% 			elseif color_scheme == 4    % ICEFERNO
% 
% 				% Make custom color schemes
% 				clrs = [1, 1, 1;
% 					0 1 1;
% 					1 0 1;
% 					1 1 0;
% 					0 1 0];
% 				nr_colors = timediff;
% 		%         nr_colors = 100;
% 				steps1 = 0.25 * nr_colors; % percent of nr of colors
% 				steps2 = 0.25 * nr_colors + steps1;
% 				steps3 = 0.25 * nr_colors + steps2;
% 				steps4 = 0.25 * nr_colors + steps3;
% 				steps = [0; steps1; steps2; steps3; steps4];
% 				background_cls = interp1(steps/nr_colors, clrs, linspace(0, 1, nr_colors));
% 		%         background_cls
% 				%background_cls = colcode.realtime(neurites.time - min(neurites.time)+1,:);
% 				line_cls = background_cls(neurites.time - min(neurites.time)+1,:);
% 			elseif color_scheme == 7	% Plain white
% 				line_cls = repmat([1,1,1], size(neurites,1), 1);
% 				background_cls = line_cls;
% 			end
			
			
			nr_colors = 50;

			if color_scheme == "Speed"  
				col_var = 'speed_ms';
				clrs = [1 0 0;  
						1 1 0;
						1 1 1];
				steps1 = 0.5 * nr_colors; % percent of nr of colors
				steps2 = 0.5 * nr_colors + steps1;
				steps = [0; steps1; steps2];
				background_cls = interp1(steps/nr_colors, clrs, linspace(0, 1, nr_colors));
			elseif color_scheme == "Amplitude"
				col_var = 'amplitudes_log';
				clrs = [1 0 0;
						1 1 0;  
						1 1 1; % I realise this is misleading but I made it like this because of the huge initial amplitude
						1 1 1];
				steps1 = 0.5 * nr_colors; % percent of nr of colors
				steps2 = 0.5 * nr_colors + steps1;
				steps3 = 0.5 * nr_colors + steps2;
				steps = [0; steps1; steps2; steps3];
				background_cls = interp1(steps/nr_colors, clrs, linspace(0, 1, nr_colors));
				neurites.amplitudes_log = log(neurites.amplitude);
			elseif  color_scheme == "Time 1" % Time
				col_var = 'time';
				% ICEFERNO colorscheme
				clrs = [1 1 1;  
						0 1 1;
						1 0 1;
						1 1 0;
						0 1 0];
				steps1 = 0.25 * nr_colors; % percent of nr of colors
				steps2 = 0.25 * nr_colors + steps1;
				steps3 = 0.25 * nr_colors + steps2;
				steps4 = 0.25 * nr_colors + steps3;
				steps = [0; steps1; steps2; steps3; steps4];
				background_cls = interp1(steps/nr_colors, clrs, linspace(0, 1, nr_colors));
			else % Time (with more colors)
				col_var = 'time';
				clrs = [0.5 0.5 0.5; 
						1 1 1;
						0 0 1;
						0 1 1;
						1 0 1;
						1 1 0;
						0 1 0;
						1 0 0];
				stepsize = 1/(size(clrs,1)-1);
				steps1 = stepsize * nr_colors; % percent of nr of colors
				steps2 = stepsize * nr_colors + steps1;
				steps3 = stepsize * nr_colors + steps2;
				steps4 = stepsize * nr_colors + steps3;
				steps5 = stepsize * nr_colors + steps4;
				steps6 = stepsize * nr_colors + steps5;
				steps = [0; steps1; steps2; steps3; steps4; steps5; steps6; nr_colors];
				background_cls = interp1(steps/nr_colors, clrs, linspace(0, 1, nr_colors));
			end
			%nr_colors = max(app.neurites{:,col_var} - min(app.neurites{:,col_var})+1);
			col_var_min = min(neurites{:,col_var});
			col_var_max = max(neurites{:,col_var});

			try % If for some reason we cant calculate the colors correctly
				line_cls = background_cls(round(((neurites{:,col_var}-col_var_min) * (1/(col_var_max - col_var_min))) * (nr_colors-1)) +1,:);
			catch
				line_cls = repmat([0.6,0.6,0.6], size(neurites,1), 1);
			end
			

			time_range = [1:size(dt_filtered,3)];% max((obj.GetSignalInit().time-5),1):(min(max(peaks.time) + 5, size(obj.GetDtTracking(),3))); % 5 time frames before/after the first/last peak;
			if ~isempty(alt_time_range)
				time_range = alt_time_range;
			end
			
			movie_frames = repmat(struct('cdata',1,'colormap',2), size(time_range,2), 1 );	% I think there's a memory leak in matlabs writeVideo function. See if it helps if we pass a struct with all the frames instead of individual frames
			i_mf = 1;

			for i_tr = time_range

				myFrame = 1;
				mytime = myFrame/10;
				fgr = figure('name','Action-Potential Propagation','NumberTitle','off', 'color', [1 1 1]*.3, 'visible','off');
				colcode.realtime = gray;
				subplot(1,1,1);
				xlim([0 electrodes_x]);
				ylim([0 electrodes_y]);
				daspect([1 1 1]);
				axis xy;
				figsizex = 1500;
				figsizey = 800;
				figposx = 500;
				figposy = 400;
				set(gcf,'Position',[figposx figposy figsizex figsizey]);
				set(gca, 'FontSize', 8)
				
				dt_filtered_i = dt_filtered(:,:,i_tr);

				% logscale
				dt_filtered_i = dt_filtered_i+0.1;
				dt_filtered_i(dt_filtered_i <= 0) = 0.00001;
				dt_filtered_i_log = log(dt_filtered_i);

				mytime = i_tr/10;

				% Show the heatmap
% 				cla(fgr);
				hold on; 

				imagesc(1:electrodes_x, 1:electrodes_y, dt_filtered_i_log);
				text(1,electrodes_y+2, num2str(i_tr));

				% Show the initiation sites
				if i_tr >= signal_init.time
					if color_scheme ~= "None" & ~plot_special_peaks
% 						scatter(mypanel.sub3(1), signal_init.x, signal_init.y, 150, 'filled','white');%,cls,'filled','s'); 
						scatter(signal_init.x, signal_init.y, 150, 'o','white','LineWidth',3);%,cls,'filled','s');
					end
					%scatter(mypanel.sub3(1), init_extras.other_initiation_sites_y, init_extras.other_initiation_sites_x, 50, 'filled','s'); 
					%hold on;
					%scatter(init_extras.initiation_sites_y, init_extras.initiation_sites_x, 1, 'filled','s'); 
					%hold on;
				end

				% Plot lines
				if i_tr > signal_init.time
					if ~isempty(neurites)
						neurites_i = neurites(neurites.time <= i_tr,:);

						if color_scheme ~= "None"
							% Lines plot neurites
							lines_plot = line([neurites_i.x_from'; neurites_i.x_to'], [neurites_i.y_from'; neurites_i.y_to']);    %,'colororder',line_cls);

							% A hack needed to differently color individual lines
							if plot_special_peaks
								line_cls(:,4) = 1; % Add Alpha
							end
							for k = 1 : length(lines_plot)
								set(lines_plot(k), {'Color','LineWidth'}, {line_cls(k,:), 1.2});
							end
						end
						
						if plot_special_peaks
							peaks_up_till_now = peaks(peaks.time <= i_tr,:);
							scatter_colors = repmat([1,1,1], size(peaks_up_till_now,1),1);
							end_point_idx = ismember(peaks_up_till_now{:,["x", "y", "time"]}, neurites{neurites.is_end_point == 1, ["x_to", "y_to", "time"]}, 'rows');
							branch_point_idx = ismember(peaks_up_till_now{:,["x", "y", "time"]}, neurites{neurites.will_branch == 1, ["x_to", "y_to", "time"]}, 'rows');
							init_idx = ismember(peaks_up_till_now{:,["x", "y", "time"]}, [signal_init.x, signal_init.y, signal_init.time], 'rows');
							new_point_idx = ismember(peaks_up_till_now.id, neurites{~ismember(neurites{:,'peak_id_from'}, neurites{:,'peak_id_to'}), 'peak_id_from'},'rows');
							new_point_idx = (new_point_idx - init_idx)==1; % Exclude the initiation site (==1 -> convert to boolean)
							normal_idx = ~(end_point_idx | branch_point_idx | init_idx);
							special_idx = (end_point_idx | branch_point_idx | init_idx | new_point_idx);
							scatter_colors(end_point_idx,:) = repmat([0.2,0.7,0.2], sum(end_point_idx), 1);
							scatter_colors(branch_point_idx,:) = repmat([1,1,1], sum(branch_point_idx), 1);
							scatter_colors(init_idx,:) = repmat([1,1,0.2], sum(init_idx), 1);
							scatter_colors(new_point_idx,:) = repmat([1,0.2,0.2], sum(new_point_idx), 1);

							scatter_size = repmat(40, size(peaks_up_till_now,1), 1);
							scatter_size(init_idx,:) = repmat(80, sum(init_idx), 1);

							% Special peaks (endpoints, branching points, initiation site, new points)
							scatter(peaks_up_till_now.x(special_idx), peaks_up_till_now.y(special_idx), scatter_size(special_idx,:), scatter_colors(special_idx,:), 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
						end
						if plot_peaks
							% Scatter plot peaks
% 							a = scatter(peaks{peaks.time <= i_tr,'x'}, peaks{peaks.time <= i_tr,'y'}, 12, 'white', 'LineWidth',0.2, 'MarkerEdgeAlpha', 0.6);
							a = scatter(peaks{peaks.time <= i_tr,'x'}, peaks{peaks.time <= i_tr,'y'}, 8, 'filled', 'white');
						end
						


					end
				end


				if plot_special_peaks
					% Legend for the special peaks
					% Need a hacky way because matlab is #$^$%&!#$^&@%#%$&#@&!
					h = zeros(3, 1);
					h(1) = scatter(NaN,NaN, 'filled', 'MarkerFaceColor', [1,1,0.2], 'DisplayName','Initiation site');
					h(2) = scatter(NaN,NaN, 'filled', 'MarkerFaceColor', [1,1,1], 'DisplayName','Branching point');
					h(3) = scatter(NaN,NaN, 'filled', 'MarkerFaceColor', [0.2,0.7,0.2], 'DisplayName','End point');
					h(4) = scatter(NaN,NaN, 'MarkerEdgeColor', [1,1,1], 'DisplayName','Data point');
					legend(h,'Location','northwest');
					
% 					c = colorbar;
% 					colormap(background_cls);
% 					c.Limits = ([0, timediff]*50);
% 					caxis([0, timediff]*50);
% 					c.Location = 'manual';
% 					c.Position = [0.1, 0.059, 0.01, 0.3];
% 					c.FontSize = 6;
% 					c.Title.String = "time (µs)";
% 					c.Title.Color = [1,1,1];
% 					c.Title.HorizontalAlignment = 'left';
				end


				%text(mypanel.sub3(1), 10, 2120, ['time: '  num2str(mytime, '%.4f') ' ms']    ,'fontsize', 9, 'HorizontalAlignment', 'left', 'color', 'k');
				colormap(inferno()); 
				caxis([-2.5 1.2]); 
				drawnow;

				movie_frames(i_mf) = getframe(fgr);
				i_mf = i_mf + 1;

				% Write gif
				if make_gif
					set(gca,'box','off')
					axis off;
					F2 = getframe(fgr);
					im = frame2im(getframe(fgr)); 
					[imind,cm] = rgb2ind(im,256); 
					% Write to the GIF File 
					if i_tr == time_range(1) 
						imwrite(imind,cm,gifName,'gif', 'Loopcount',Inf, 'DelayTime', 0.1); 
					elseif i_tr == time_range(end)
						imwrite(imind,cm,gifName,'gif','WriteMode','append','DelayTime',3); 
						imwrite(imind,cm,strcat(gifName,'.png'),'png'); 
					else
						imwrite(imind,cm,gifName,'gif','WriteMode','append','DelayTime',0.1); 
					end
					set(gca,'box','on')
					axis on;
				end

				close(fgr);
			end
			
			
				
			% Write video
			% Because somehow I randomly get this error: "Frame must be 1500 by 800". IT IS 1500 x 800 !!!
			myVideo = VideoWriter(videoName); % If you get an error here you might need to run: close(myVideo)
			myVideo.FrameRate = 10;
			myVideo.Quality = 100;
			open(myVideo);
			errors = 0;
			while errors < 10 % just try 10 times and other wise move on
				try
					writeVideo(myVideo,movie_frames);
					break;
				catch E
					errors = errors + 1;
				end
			end
			close(myVideo);
			
			
			

			if make_eps
				%plot2svg(svgName, lines_plot)

				make_it_tight = true;
				subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0],0.06,0.031); % No borders/margins: subplot = @(m,n,p) subtightplot (m, n, p, [0.00 0.00],0,0);
				if ~make_it_tight,  clear subplot;  end

				f = figure('color','black', 'Visible','off');
				subplot(1,1,1);
				
				if plot_special_peaks == 1
					scatter_colors = repmat([1,1,1], size(peaks,1),1);
					end_point_idx = ismember(peaks{:,["x", "y", "time"]}, neurites{neurites.is_end_point == 1, ["x_to", "y_to", "time"]}, 'rows');
					branch_point_idx = ismember(peaks{:,["x", "y", "time"]}, neurites{neurites.will_branch == 1, ["x_to", "y_to", "time"]}, 'rows');
					init_idx = ismember(peaks{:,["x", "y", "time"]}, [signal_init.x, signal_init.y, signal_init.time], 'rows');
					new_point_idx = ismember(peaks.id, neurites{~ismember(neurites{:,'peak_id_from'}, neurites{:,'peak_id_to'}), 'peak_id_from'},'rows');
					new_point_idx = (new_point_idx - init_idx)==1; % Exclude the initiation site (==1 -> convert to boolean)
					normal_idx = ~(end_point_idx | branch_point_idx | init_idx);
					special_idx = (end_point_idx | branch_point_idx | init_idx | new_point_idx);
					scatter_colors(end_point_idx,:) = repmat([0.2,0.7,0.2], sum(end_point_idx), 1);
					scatter_colors(branch_point_idx,:) = repmat([1,1,1], sum(branch_point_idx), 1);
					scatter_colors(init_idx,:) = repmat([1,1,0.2], sum(init_idx), 1);
					scatter_colors(new_point_idx,:) = repmat([1,0.2,0.2], sum(new_point_idx), 1);

					scatter_size = repmat(40, size(peaks,1), 1);
					scatter_size(init_idx,:) = repmat(80, sum(init_idx), 1);

					% Special peaks (endpoints, branching points, initiation site, new points)
					scatter(peaks.x(special_idx), peaks.y(special_idx), scatter_size(special_idx,:), scatter_colors(special_idx,:), 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
					hold on;
				else
					% No special peaks, but at least show the initiation site
					scatter(signal_init.x, signal_init.y, 150, 'filled', 'o','white','LineWidth',3);%,cls,'filled','s'); 
					hold on;
				end
				if plot_peaks == 1
					scatter(peaks.x, peaks.y, 6, 'filled', 's', 'white');
				end
				
				
				if color_scheme ~= "None"
					lines_plot = line([neurites.x_from'; neurites.x_to'], [neurites.y_from'; neurites.y_to']);
					hold on; 
			%         set(findobj(gcf, 'type','axes'), 'Visible','off');
					% A hack needed to differently color individual lines
					line_cls(:,4) = 1; % Add Alpha
					for k = 1 : length(lines_plot)
						set(lines_plot(k), {'Color','LineWidth'}, {line_cls(k,:), 3});
					end

					c = colorbar;
					colormap(background_cls);
					c.Limits = ([0, timediff]*50);
					caxis([0, timediff]*50);
					c.Location = 'east';
					c.Position = [0.95, 0.059, 0.01, 0.3];
					c.FontSize = 6;
					c.Title.String = "time (µs)";
					c.Title.Color = [1,1,1];
					c.Title.HorizontalAlignment = 'left';
				end

				if plot_special_peaks == 1
					% Legend for the special peaks
					% Need a hacky way because matlab is #$^$%&!#$^&@%#%$&#@&!
					h = zeros(3, 1);
					hold on;
					h(1) = scatter(NaN,NaN, 'filled', 'MarkerFaceColor', [1,1,0.2], 'DisplayName','Initiation site');
					h(2) = scatter(NaN,NaN, 'filled', 'MarkerFaceColor', [1,1,1], 'DisplayName','Branching point');
					h(3) = scatter(NaN,NaN, 'filled', 'MarkerFaceColor', [0.2,0.7,0.2], 'DisplayName','End point');
					h(4) = scatter(NaN,NaN, 'MarkerEdgeColor', [1,1,1], 'DisplayName','Data point');
					legend(h);
				end

		        whitebg('black');
				xlim([0,220]);
				ylim([0,120]);
				xlabel('X-coordinates (electrodes)');
				ylabel('Y-coordinates (electrodes)');
				set(gcf, 'color', [0,0,0]);
				daspect([1 1 1]);

				f.InvertHardcopy = 'off'; % Preserve colors when saving
				rez=200; %resolution (dpi) of final graphic
				%f=f1; %f is the handle of the figure you want to export
				figpos=getpixelposition(f); %dont need to change anything here
				resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
				set(f,'paperunits','points','papersize',figpos(3:4)/resolution,...
				'paperposition',[0 0 220*4 120*4]); %paperposition can be used to actually increase figure size

				print(f,pngName,'-dpng',['-r',num2str(rez)],'-opengl'); %save file
				print(f,svgName, '-depsc');
				close(f);
			end
			
			%disp(strcat('output/movies/',dirname,'_',dfile,' - DONE'));
		end
		
		function nul = ExportData(obj)
			obj.Output("ExportData()", 0);
			
			file_base_name = strcat(obj.GetRecordingId().path_tracking);

			neurites = obj.GetNeurites();
			peaks = obj.GetPeaks();
			settings = obj.GetSettings();
			dt_filtered = obj.GetDtTracking();
			
			save(strcat(file_base_name,'.neurites.mat'), 'neurites');
			save(strcat(file_base_name,'.settings.mat'), 'settings');
			save(strcat(file_base_name,'.peaks.mat'), 'peaks');
			save(strcat(file_base_name,'.dt_filtered.mat'), 'dt_filtered');
		end
		
		function ret = CreateEmptyTableCopy(obj, original_table, n_rows)
			% Make a local empty copy of a table
			empty_table_copy = original_table;	% COPY
			empty_table_copy(:,:) = [];	% EMPTY (0 rows)
			empty_table_copy{1:n_rows,:} = NaN; % Create rows with NaN values
			ret = empty_table_copy;
		end
		
		function nul = ConnectPeaksManually(obj, peak_id_a, peak_id_b)
			obj.Output("ConnectPeaksManually()", 0);
			
			peaks = obj.GetPeaks();
			peak_a = peaks(peaks.id == peak_id_a,:);
			peak_b = peaks(peaks.id == peak_id_b,:);
			
			if peak_a.time == peak_b.time
				disp("Warning - Not connecting peaks, they occur at the same time");
				return
			end
			
			if peak_a.time < peak_b.time
				peak_from = peak_a;
				peak_to = peak_b;
			else
				peak_from = peak_b;
				peak_to = peak_a;
			end
			
			% Distance in electrodes
			distance_el = abs(peak_from{1,{'x','y'}} - peak_to{1,{'x','y'}});
			distance_el = sqrt(distance_el(:,1).^2 + distance_el(:,2).^2);
			
			new_neurites = obj.CreateEmptyTableCopy(obj.GetNeurites(),1);
			new_neurites.x_from = peak_from.x;
			new_neurites.y_from = peak_from.y;
			new_neurites.x_to = peak_to.x;
			new_neurites.y_to = peak_to.y;
			new_neurites.time = peak_to.time;
			new_neurites.speed_el_tf = distance_el / (peak_to.time - peak_from.time); % Speed in electrodes per timeframe
			new_neurites.length_el = distance_el;
			new_neurites.speed_ms = new_neurites.speed_el_tf * (17.5/50); % pitch (17.5 um) / duration of time frame (50 us)
			new_neurites.amplitude = peak_to.amplitude;
			new_neurites.voltage = peak_to.voltage;
			new_neurites.peak_id_from = peak_from.id;
			new_neurites.peak_id_to = peak_to.id;
			new_neurites.track_type = 4; % 4 = manual connection
			new_neurites.previous_peak_time = peak_from.time;
			
			obj.AddNeurites(new_neurites);
			
			obj.AppendNeuriteLengthsAndIdsAndStuff();
		end
		
		function nul = LoadTrack(obj, recording_id)
			obj.Output("LoadTrack()", 0);
			
			obj.SetRecordingId(recording_id);
			obj.ClearVariables();
			obj.LoadExportedData();
		end
		
		function nul = NewTrack(obj, recording_id)
			obj.Output("NewTrack()", 0);
			
			obj.SetRecordingId(recording_id);
			obj.ClearVariables();
		end
		
		function nul = ReconnectPeaks(obj)
			obj.Output("ReconnectPeaks()", 0);
			
			obj.neurites = [];
			
			obj.ConnectPeaks();
		end
		
		function ret = ReloadDtAmplitude(obj)
			% Recreate dt_amplitude, in case we changed some settings
			obj.dt_amplitude = [];
			ret = obj.GetDtAmplitude();
		end
		
		function ret = ReloadDtTracking(obj)
			% Recreate dt_tracking, in case we changed some settings
			obj.dt_tracking = [];
			ret = obj.GetDtTracking();
		end
		
		function nul = RemoveStationaryPeaks(obj)
			% In an attempt to remove unwanted peaks, we will filter out those that do not move.
			% These for example are found when multiple neurons spike simultaneously and their initiation sites are recorded
			
% 			peaks = obj.GetPeaks();
% 			peak_ids_to_remove = [];
% 			for i_p = 1:size(peaks,1)
% 				
% 				% if the time between the peaks is between 4 and 10, this is likely a pre or post peak
% 				peaks_close_time_range = peaks(abs(peaks.time - peaks.time(i_p)) <= 2 & peaks.id ~= peaks.id(i_p),:);
% 				
% 				if isempty(peaks_close_time_range); continue; end
% 				
% 				% Get closest peak
% 				distances =  sqrt(sum((peaks_close_time_range{:,{'x','y'}} - peaks{i_p,{'x','y'}}).^2,2));
% 				
% 				if isempty(distances(distances <= 0.5)); continue; end
% 				
% % 				close_peaks = peaks_close_time_range(distances < 0.5:);
% 				
% 				% Remove both peaks, as they appear stationary
% % 				if peaks.amplitude(i_p) < max(close_peaks.amplitude)
% 					% Delete
% 				peak_ids_to_remove = [peak_ids_to_remove; peaks.id(i_p)];
% % 				end
% 			end

			min_speed = 0.13; % meters / second

			neurites = obj.GetNeurites();

			no_inf_loop = 100;

			% Remove near stationary end points
			peak_ids_to_remove = neurites{neurites.is_end_point == 1 & neurites.speed_ms < min_speed,'peak_id_to'};
			while ~isempty(peak_ids_to_remove) && no_inf_loop > 0
				obj.DelPeaks(peak_ids_to_remove);
				obj.AppendNeuriteLengthsAndIdsAndStuff();
				
				neurites = obj.GetNeurites();	% Get updated neurites table
				peak_ids_to_remove = neurites{neurites.is_end_point == 1 & neurites.speed_ms < min_speed,'peak_id_to'};
				
				no_inf_loop = no_inf_loop - 1;
			end

			% Remove near stationary start points
			peaks = obj.GetPeaks();
			signal_init_peak_id = peaks.id(peaks.is_signal_init == 1);
			no_inf_loop = 100;
			peak_ids_to_remove = neurites{~ismember(neurites{:,'peak_id_from'}, neurites{:,'peak_id_to'}) & neurites.speed_ms < min_speed, 'peak_id_from'};
			while ~isempty(peak_ids_to_remove) && no_inf_loop > 0
				obj.DelPeaks(peak_ids_to_remove);
				obj.AppendNeuriteLengthsAndIdsAndStuff();
				
				neurites = obj.GetNeurites();	% Get updated neurites table
				peak_ids_to_remove = neurites{~ismember(neurites{:,'peak_id_from'}, neurites{:,'peak_id_to'}) & neurites.speed_ms < min_speed, 'peak_id_from'};
				peak_ids_to_remove = peak_ids_to_remove(peak_ids_to_remove ~= signal_init_peak_id);
				
				no_inf_loop = no_inf_loop - 1;
			end
			
		end
		
		function nul = RemoveTinyBranches(obj)
			% Attempt to remove "tiny" branches, these can be an artifact of some noise filtering, or can occur when the signal
			% is a bit "jagged" or jumpy

			min_branch_length_el = 2.5; % electrodes

			neurites = obj.GetNeurites();
			
			% Calculate the length of the last branch (last segment of the last branch has is_end_point = 1)
			% calculate size of that last branch by starting to calculate the length from the first segment, so right after the branching point
			last_branches = neurites(ismember(neurites(:,"neurite_id"), neurites(neurites.is_end_point == 1,"neurite_id")) & neurites.segment_id == 1, :);
			last_branches.last_branch_length = last_branches.branch_length - last_branches.cum_length + last_branches.length_el;
			too_short_last_branches = last_branches(last_branches.last_branch_length < min_branch_length_el,:);
			
			obj.DelPeaks(neurites{ismember(neurites.neurite_id, too_short_last_branches.neurite_id), 'peak_id_to'});
		end
	end

	%% ----- METHODS (PRIVATE) -----
	methods (Access = private)

		function nul = AppendNeuriteLengthsAndIdsAndStuff(obj)
			obj.Output("AppendNeuriteLengthsAndIdsAndStuff()", 0);
			
			% Start at the beginning: neuriteID = 1; Level = 1; SegmentID = 1
			% Find all connections:
			% If 1 connection
			% NeuriteID = same; Level = same; SegmentID = parent + 1;
			% If >1 connections
			% NeuriteID = Parent + [1,2,3]; Level = parent + 1; SegmentID = 1;
			
			%neurites2 = dataset({pc_neurites, 'x_from','y_from','x_to','y_to','time','line_id','speed','length','speed_ms','amplitude', 'peak_id_from', 'peak_id_to'});
			neurites2 = obj.GetNeurites();
			neurites2.neurite_id(:) = "-";
			neurites2.segment_id(:) = 0;
			neurites2.branch_level(:) = 0;     % Number of branch points we passed before getting to this segment
			neurites2.cum_length(:) = 0;         % Total distance from the initiation site to the end of that segment
			neurites2.is_end_point(:) = 0;     % Boolean: Segment does not connect to other points Y/N?
			neurites2.will_branch(:) = 0;      % Boolean: End of segment branches Y/N?
			neurites2.branch_length(:) = 0;    % Length from the initiation site to that segment, PLUS all other segments BEYOND that point

			% Get level 1 segments
			l1_seg_idx = find(neurites2.x_from == obj.GetSignalInit().x & neurites2.y_from == obj.GetSignalInit().y);

			neurite_id_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

			new_neurite_id = 1;
% 			disp(l1_seg_idx)
			for i = l1_seg_idx'
				neurites2.neurite_id(i) = convertCharsToStrings(neurite_id_chars(new_neurite_id));
				neurites2.segment_id(i) = 1;
				neurites2.branch_level(i) = 0;
				neurites2.cum_length(i) = neurites2.length_el(i);
				new_neurite_id = new_neurite_id + 1;
				
			end

			% Now loop over each main branch?
			parent_ids = find(neurites2.neurite_id ~= "-");
			we_have_daughters = 1;
			while we_have_daughters
				new_parent_ids = [];
				for i = parent_ids'
					parent_segment = neurites2(i,:);
					daughter_ids = find(neurites2.x_from == parent_segment.x_to & neurites2.y_from == parent_segment.y_to & neurites2.neurite_id == "-" & neurites2.previous_peak_time == parent_segment.time);
					if size(daughter_ids,1) == 1
						neurites2.neurite_id(daughter_ids) = parent_segment.neurite_id;
						neurites2.segment_id(daughter_ids) = parent_segment.segment_id + 1;
						neurites2.branch_level(daughter_ids) = parent_segment.branch_level;
						neurites2.cum_length(daughter_ids) = parent_segment.cum_length + neurites2.length_el(daughter_ids);
					elseif size(daughter_ids,1) > 1
						new_neurite_id = 1;
						for j = daughter_ids'
							%disp('a')
							%parent_segment.neurite_id
							%new_neurite_id
							new_neurite_id2 = strcat(parent_segment.neurite_id, neurite_id_chars(new_neurite_id));
							neurites2.neurite_id(j) = convertCharsToStrings(new_neurite_id2);
							neurites2.segment_id(j) = 1;
							neurites2.branch_level(j) = parent_segment.branch_level+1;
							neurites2.cum_length(j) = parent_segment.cum_length + neurites2.length_el(j);
							new_neurite_id = new_neurite_id + 1;
							%disp(strcat(num2str(parent_segment.neurite_id, '%.0f'), " - ", num2str(new_neurite_id, '%.0f')));
						end
						neurites2.will_branch(i) = 1;
					elseif size(daughter_ids,1) == 0
						neurites2.is_end_point(i) = 1;
					end
					new_parent_ids = [new_parent_ids; daughter_ids];
				end

				parent_ids = new_parent_ids;

				if isempty(parent_ids)
					we_have_daughters = 0;
				end
			end

		%     % Don't need this anymore since neurite_id is a string now. Keeping it for future reference
		%     max_branch_level = max(neurites2.branch_level);
		%     for i = 1:size(neurites2,1)
		%         n = neurites2.neurite_id(i);
		%         neurites2.neurite_id(i) = (n * 10^((max_branch_level+1) - size(num2str(n, '%.0f'),2)));
		%     end


			max_branch_level = max(neurites2.branch_level);
			for i = 1:size(neurites2,1)
				neurite = neurites2(i,:);
				if ismissing(neurite.neurite_id)
					continue;
				end
				
				lower_segments = neurites2((startsWith(neurites2.neurite_id, neurite.neurite_id) & neurites2.neurite_id ~= neurite.neurite_id) | (neurites2.neurite_id == neurite.neurite_id & neurites2.segment_id > neurite.segment_id), :);
				
	%         % Version I used for the numerical neurite_id:
	%         lower_segments = neurites2((neurites2.neurite_id > neurite.neurite_id & neurites2.neurite_id < (neurite.neurite_id + 1 * 10^(max_branch_level - neurite.branch_level))...
	%             | (neurites2.neurite_id == neurite.neurite_id & neurites2.segment_id > neurite.segment_id)),:);
				lower_segments_length_sum = sum(lower_segments.length_el);


				neurites2.branch_length(i) = neurite.cum_length + lower_segments_length_sum;
			end
			
			obj.neurites = neurites2;
		end
		
		function nul = CalculateAbsoluteNoiseThresholds(obj)
			settings = obj.GetSettings();
			dt_tracking = obj.GetDtTracking();
			
			if settings.noise_threshold.use_plain_thresh_test == 1 || obj.GetNeuronType() == obj.neuron_types(3)
				obj.Output("Neuron type unknown, using plain noise threshold settings", 0);
				threshold_std = settings.noise_threshold.plain_thresh_test.std;
			else
				threshold_std = settings.noise_threshold.(obj.GetNeuronType()).std;
			end
			
			settings.noise_threshold.(obj.GetNeuronType()).absolute.high = std(dt_tracking(:,:,1:10), 0, 'all') * threshold_std.high * settings.noise_threshold_modifier;
			settings.noise_threshold.(obj.GetNeuronType()).absolute.medium = std(dt_tracking(:,:,1:10), 0, 'all') * threshold_std.medium * settings.noise_threshold_modifier;
			settings.noise_threshold.(obj.GetNeuronType()).absolute.low = std(dt_tracking(:,:,1:10), 0, 'all') * threshold_std.low * settings.noise_threshold_modifier;
			
			obj.settings = settings;
		end
		
		function ret = CalculateNeuronStats(obj)
			obj.Output("CalculateNeuronStats()", 0);
			
			neuron_stats = obj.GetNeuronStats();
			
			neurites = obj.GetNeurites();
			peaks = obj.GetPeaks();
			signal_init = obj.GetSignalInit();
			settings = obj.GetSettings();
			
            neuron_subtypes = ["motor","cortical"];
            inter_branch_lengths = grpstats(neurites, {'neurite_id'},{'sum'},'datavars','length_el');
            
            amplitude_quartiles = quantile(peaks.amplitude(peaks.time > (signal_init.time+8)), 3);
            amplitude_deciles = quantile(peaks.amplitude(peaks.time > (signal_init.time+8)), 9);
            voltage_quartiles = quantile(peaks.voltage(peaks.time > (signal_init.time+8)), 3);
            voltage_deciles = quantile(peaks.voltage(peaks.time > (signal_init.time+8)), 9);
            speed_quartiles = quantile(neurites.speed_ms, 3);
            speed_deciles = quantile(neurites.speed_ms, 9);
            length_to_endpoints_quartiles = quantile(neurites.cum_length(neurites.is_end_point==1)*17.5, 3);
            inter_branch_length_quartiles = quantile(inter_branch_lengths.sum_length_el*17.5, 3);
            time_to_endpoints_quartiles = quantile((neurites.time(neurites.is_end_point==1) - signal_init.time)*50, 3);
			
            %neuron_stats.recording_id = convertCharsToStrings(obj.GetRecordingId().recording_id);
			if ismissing(neuron_stats.track_quality)
				neuron_stats.track_quality = "unknown";
			end
			neuron_stats.axtracktor_version = convertCharsToStrings(settings.axtracktor_version);
			neuron_stats.subtype = obj.GetNeuronType(); 
            neuron_stats.amplitude_initiation = peaks.amplitude(peaks.x == signal_init.x & peaks.y == signal_init.y & peaks.time == signal_init.time);
            neuron_stats.amplitude_mean = mean(peaks.amplitude(peaks.time > (signal_init.time+8)));
            neuron_stats.amplitude_std = std(peaks.amplitude(peaks.time > (signal_init.time+8)));
            neuron_stats.amplitude_median = amplitude_quartiles(2);
            neuron_stats.amplitude_iqr = iqr(peaks.amplitude(peaks.time > (signal_init.time+8)));
            neuron_stats.amplitude_q1 = amplitude_quartiles(1);
            neuron_stats.amplitude_q3 = amplitude_quartiles(3);
            neuron_stats.amplitude_d1 = amplitude_deciles(1);
            neuron_stats.amplitude_d2 = amplitude_deciles(2);
            neuron_stats.amplitude_d3 = amplitude_deciles(3);
            neuron_stats.amplitude_d4 = amplitude_deciles(4);
            neuron_stats.amplitude_d6 = amplitude_deciles(6);
            neuron_stats.amplitude_d7 = amplitude_deciles(7);
            neuron_stats.amplitude_d8 = amplitude_deciles(8);
            neuron_stats.amplitude_d9 = amplitude_deciles(9);
            neuron_stats.voltage_initiation = peaks.voltage(peaks.x == signal_init.x & peaks.y == signal_init.y & peaks.time == signal_init.time);
            neuron_stats.voltage_mean = mean(peaks.voltage(peaks.time>(signal_init.time+8)));
            neuron_stats.voltage_std = std(peaks.voltage(peaks.time > (signal_init.time+8)));
            neuron_stats.voltage_median = voltage_quartiles(2);
            neuron_stats.voltage_iqr = iqr(peaks.voltage(peaks.time > (signal_init.time+8)));
            neuron_stats.voltage_q1 = voltage_quartiles(1);
            neuron_stats.voltage_q3 = voltage_quartiles(3);
            neuron_stats.voltage_d1 = voltage_deciles(1);
            neuron_stats.voltage_d2 = voltage_deciles(2);
            neuron_stats.voltage_d3 = voltage_deciles(3);
            neuron_stats.voltage_d4 = voltage_deciles(4);
            neuron_stats.voltage_d6 = voltage_deciles(6);
            neuron_stats.voltage_d7 = voltage_deciles(7);
            neuron_stats.voltage_d8 = voltage_deciles(8);
            neuron_stats.voltage_d9 = voltage_deciles(9);
            neuron_stats.speed_mean = mean(neurites.speed_ms);
            neuron_stats.speed_std = std(neurites.speed_ms);
            neuron_stats.speed_median = speed_quartiles(2);
            neuron_stats.speed_iqr = iqr(neurites.speed_ms);
            neuron_stats.speed_q1 = speed_quartiles(1);
            neuron_stats.speed_q3 = speed_quartiles(3);
            neuron_stats.speed_d1 = speed_deciles(1);
            neuron_stats.speed_d2 = speed_deciles(2);
            neuron_stats.speed_d3 = speed_deciles(3);
            neuron_stats.speed_d4 = speed_deciles(4);
            neuron_stats.speed_d6 = speed_deciles(6);
            neuron_stats.speed_d7 = speed_deciles(7);
            neuron_stats.speed_d8 = speed_deciles(8);
            neuron_stats.speed_d9 = speed_deciles(9);
            neuron_stats.nr_endpoints = sum(neurites.is_end_point);
            neuron_stats.nr_branching_points = sum(neurites.will_branch);
            neuron_stats.total_length = sum(neurites.length_el)*17.5;
            neuron_stats.length_to_endpoints_mean = mean(neurites.cum_length(neurites.is_end_point==1)*17.5);
            neuron_stats.length_to_endpoints_std = std(neurites.cum_length(neurites.is_end_point==1)*17.5);
            neuron_stats.length_to_endpoints_median = length_to_endpoints_quartiles(2);
            neuron_stats.length_to_endpoints_iqr = iqr(neurites.cum_length(neurites.is_end_point==1)*17.5);
            neuron_stats.length_to_endpoints_q1 = length_to_endpoints_quartiles(1);
            neuron_stats.length_to_endpoints_q3 = length_to_endpoints_quartiles(3);
            neuron_stats.inter_branch_length_mean = mean(inter_branch_lengths.sum_length_el*17.5);
            neuron_stats.inter_branch_length_std = std(inter_branch_lengths.sum_length_el*17.5);
            neuron_stats.inter_branch_length_median = inter_branch_length_quartiles(2);
            neuron_stats.inter_branch_length_iqr = iqr(inter_branch_lengths.sum_length_el*17.5);
            neuron_stats.inter_branch_length_q1 = inter_branch_length_quartiles(1);
            neuron_stats.inter_branch_length_q3 = inter_branch_length_quartiles(3);
            neuron_stats.time_to_endpoints_mean = mean((neurites.time(neurites.is_end_point==1) - signal_init.time)*50);
            neuron_stats.time_to_endpoints_std = std((neurites.time(neurites.is_end_point==1) - signal_init.time)*50);
            neuron_stats.time_to_endpoints_median = time_to_endpoints_quartiles(2);
            neuron_stats.time_to_endpoints_iqr = iqr((neurites.time(neurites.is_end_point==1) - signal_init.time)*50);
            neuron_stats.time_to_endpoints_q1 = time_to_endpoints_quartiles(1);
            neuron_stats.time_to_endpoints_q3 = time_to_endpoints_quartiles(3);
			
			ret = neuron_stats;
		end
		
		function nul = ClearTrackingData(obj)
			obj.Output("ClearTrackingData()", 0);
			
			obj.peaks = [];
			obj.neurites = [];
		end
		
		function nul = ConnectPeaks(obj)
			obj.Output("ConnectPeaks()", 0);
			
			time_frames = (obj.GetSignalInit().time+1):size(obj.GetDtTracking(),3); % Going backwards 
			
			for time_frame = time_frames
				obj.AddNeurites(obj.ConnectPeaksInTimeFrame(time_frame));
			end
		end
		
		function ret = ConnectPeaksInTimeFrame(obj, time_frame)    
			ret = [];
			
			Angles360 = @(a) rem(360+a, 360); % kudos to Star Strider @ https://nl.mathworks.com/matlabcentral/answers/286805-transform-angles-that-are-180-to-180-to-angles-that-are-in-the-0-360-degree-range

			peaks = obj.GetPeaks();
			if ~isempty(peaks)
				previous_peaks_overall = peaks(peaks.time<time_frame,:);
				current_peaks = peaks(peaks.time == time_frame,:);
			end

			new_neurites = obj.CreateEmptyTableCopy(obj.GetNeurites(),0);
			
			if isempty(previous_peaks_overall) % Nothing to connect to/from
				ret = new_neurites;
				return;
			end
			
			electrodes_x = size(obj.GetDtTracking(),2);
			electrodes_y = size(obj.GetDtTracking(),1);

			current_neurite_id = 1;
			for j=1:size(current_peaks,1)
				previous_peaks = peaks(peaks.time>=time_frame-obj.GetSettings().connect_look_back_time & peaks.time<time_frame,:);

				% Get distance between all other peaks
				dist = abs(previous_peaks{:,{'x','y'}} - current_peaks{j,{'x','y'}});
				dist = sqrt(dist(:,1).^2 + dist(:,2).^2);
				[min_dist, min_index] = min(dist);
				
				if isempty(min_dist) % Might occur if no previous peaks are found
					min_dist = -1;
				end
				
				% Connect it to the closest peak
				if min_dist>0 && min_dist < obj.GetSettings().connect_direct_min_peak_dist
					previous_peak.x = previous_peaks{min_index,'x'};  % x-pos
					previous_peak.y = previous_peaks{min_index,'y'};  % y-pos
					previous_peak.time = previous_peaks{min_index,'time'};  % time frame

					i_nn = size(new_neurites,1)+1;
					new_neurites{i_nn,:} = NaN; % Initiate a new ("empty") row
					new_neurites.x_from(i_nn) = previous_peaks{min_index,'x'};
					new_neurites.y_from(i_nn) = previous_peaks{min_index,'y'};
					new_neurites.x_to(i_nn) = current_peaks{j,'x'};
					new_neurites.y_to(i_nn) = current_peaks{j,'y'};
					new_neurites.time(i_nn) = time_frame;
					new_neurites.speed_el_tf(i_nn) = min_dist / (time_frame - previous_peak.time); % Speed in electrodes per timeframe
					new_neurites.length_el(i_nn) = min_dist;
					new_neurites.speed_ms(i_nn) = new_neurites.speed_el_tf(i_nn) * (17.5/50); % pitch (17.5 um) / duration of time frame (50 us)
					new_neurites.amplitude(i_nn) = current_peaks{j,'amplitude'};
					new_neurites.voltage(i_nn) = current_peaks{j,'voltage'};
					new_neurites.peak_id_from(i_nn) = previous_peaks{min_index,'id'};
					new_neurites.peak_id_to(i_nn) = current_peaks{j,'id'}; % old: max(previous_peaks(:,'id')) + j;
					new_neurites.track_type(i_nn) = 1; % 1 = direct connection
					new_neurites.previous_peak_time(i_nn) = previous_peaks{min_index,'time'};
					
				% Peaks are too far for a direct connection, use dynamic programming
				elseif obj.GetSettings().use_dynamic_programming
					dt_tracking = obj.GetDtTracking();
					dt_voltage = obj.GetDtVoltage();
					dt_amplitude = obj.GetDtAmplitude();
					neurites = obj.GetNeurites();
					
					% Find all previous peaks within the time range
					previous_peaks = peaks(peaks.time>=time_frame-obj.GetSettings().dp.time_range & peaks.time<time_frame,:);
					% Find all previous peaks within search radius
					peak_distances = abs(previous_peaks{:,{'x','y'}} - current_peaks{j,{'x','y'}});
					peak_distances = sqrt(peak_distances(:,1).^2 + peak_distances(:,2).^2); % Simple euclidian distance
					previous_peaks = previous_peaks(peak_distances < obj.GetSettings().dp.radius,:);
					
					routes = struct();
					% For each peak in radius & time range, get some statistics
					for i_pp = 1:size(previous_peaks,1)
						% Get all electrodes in a straight line between current peak and the previous peak
						% Basically we find each electrod between these two points, and get the signal strength for those locations
						spatiotemporal_diff = current_peaks{j, {'x','y','time'}} - previous_peaks{i_pp, {'x','y','time'}};
						steps = max(abs(round(spatiotemporal_diff)));
						route_signal_strength = zeros([1,steps]);
						for i_step = 1:steps
							step_coords = round(previous_peaks{i_pp, {'x','y','time'}} + ((spatiotemporal_diff./steps)*i_step));
							route_signal_strength(i_step) = dt_tracking(step_coords(2), step_coords(1), step_coords(3)); % Note: first index is Y, second is X, (third is time)
						end
						
						routes(i_pp).signal_mean = mean(route_signal_strength);
						routes(i_pp).signal_std = std(route_signal_strength);
						tmp_signal_quantile = quantile(route_signal_strength, 3);
						routes(i_pp).signal_q1 = tmp_signal_quantile(1);
						routes(i_pp).distance = sqrt(spatiotemporal_diff(1)^2 + spatiotemporal_diff(2)^2);
						
						% Speed comparison
						previous_neurite = neurites(neurites.peak_id_to == previous_peaks.id(i_pp),:);
						if isempty(previous_neurite)
							previous_speed = 0;
						else
							previous_speed = previous_neurite.speed_el_tf; % speed in electrode per time frame
						end
						current_speed = routes(i_pp).distance / abs(spatiotemporal_diff(3));
						routes(i_pp).speed = current_speed;
						routes(i_pp).speed_diff = abs(previous_speed - current_speed);
						
						% Calculate the angle
						if isempty(previous_neurite)
							routes(i_pp).ba_angle = [];
						else
							% ab = last segment before gap
							% cd = first segment after gap
							% bc = gap (i.e. this connection we're trying to establish)
							a = [previous_neurite.x_from(1), previous_neurite.y_from(1)]; % Start point of last segment before gap
							b = [previous_neurite.x_to(end), previous_neurite.y_to(end)]; % End point of last segment before gap / Also starting point of the gap
							ba_angle = Angles360( atan2d(b(1)-a(1), b(2)-a(2)) );
							routes(i_pp).ba_angle = ba_angle;
						end
						if isempty(routes(i_pp).ba_angle)
							% Connected to starting point, so no angle can be
							% calculated. Set to 90, as an intermediate angle
							routes(i_pp).total_angle = 90;
						else
							% angle from the start to the end point
							% atan2d gets the angle using a downward vector as base line
							c = [current_peaks.x(j), current_peaks.y(j)]; % Start point of first segment after gap / Also last point of the gap
% 							d = [neurite_segment_after_gap.x_to(1), neurite_segment_after_gap.y_to(1)];
							
							% ABC angle
							bc_angle = Angles360( atan2d(b(1)-c(1), b(2)-c(2)) );
							routes(i_pp).bc_angle = bc_angle;
							angle_abc = abs( ba_angle-bc_angle );
							if angle_abc > 180
								angle_abc = 360 - angle_abc;
							end
							angle_abc_straight = 180 - angle_abc;
							 
% 							% BCD angle - In case we know what the angle of the next neurite will be (might use this in the future to improve angle calculation)
% 							cb_angle = Angles360( atan2d(c(1)-b(1), c(2)-b(2)) );
% 							cd_angle = Angles360( atan2d(c(1)-d(1), c(2)-d(2)) );
% 							routes(irpp).cb_angle = cb_angle;
% 							routes(irpp).cd_angle = cd_angle;
% 							angle_bcd = abs( cb_angle-cd_angle );
% 							if angle_bcd > 180
% 								angle_bcd = 360 - angle_bcd;
% 							end
% 							angle_bcd_straight = 180 - angle_bcd
% 							
% 							routes(irpp).bcd_angle = angle_abc_straight; % 0 is straight ahead, 180 is going back in the direction we came
							
							% Total angle (abc + bcd)
							routes(i_pp).total_angle = angle_abc_straight; % + angle_bcd_straight
						end
					end
					
					if ~isfield(routes, 'signal_mean') || size(routes,2) == 0 || isempty(routes)	% isempty() doesn't work on "empty" structs, just check if one of the fieldnames exists or not instead
						continue; % No connection possible
					end
					
					% Try some combined scoring metric, basically score the best on
					% average over these properties. Maybe I should make it
					% gaussian? Also they could be weighted differently of course
					rank_mean = [routes.signal_mean] / max([routes.signal_mean]);	% Normalise
					quant_1 = [routes.signal_q1] + (min(0, min([routes.signal_q1]))*-1); % if we have negative values, increase them to all be 0 or positive, otherwise just leave it as is
					rank_quant_1 = repelem(1, size(rank_mean,2));
					if max(quant_1) ~= 0 % Can be 0 if they all have the same value
						rank_quant_1 = quant_1 / max(quant_1);		% Normalise
					end
					rank_std = 1 - [routes.signal_std] / max([routes.signal_std]);				% Normalise
					rank_speed_diff = 1 - [routes.speed_diff] / max([routes.speed_diff]);		% Normalise
					rank_distance = 1 - [routes.distance] / max([routes.distance]);				% Normalise
					rank_angle = 1 - [routes.total_angle] / 180;								% Not a normlisation but we want a value between 0 and 1

					% Assign weights to the different statistics and select the peak with the highest score as our connected peak
					rw = obj.GetSettings().dp.rank_weights;            
					rank = (rank_mean*rw.mean) + (rank_quant_1*rw.q_1) + (rank_std*rw.std) + (rank_speed_diff*rw.speed_diff) + (rank_distance*rw.distance) + (rank_angle*rw.angle); % Multiplications are used to increase the weight of the feature (higher is more important)
					[best_score, best_idx] = max(rank);

					i_nn = size(new_neurites,1)+1;
					new_neurites{i_nn,:} = NaN;
					new_neurites.x_from(i_nn) = previous_peaks.x(best_idx);
					new_neurites.y_from(i_nn) = previous_peaks.y(best_idx);
					new_neurites.x_to(i_nn) = current_peaks.x(j);
					new_neurites.y_to(i_nn) = current_peaks.y(j);
					new_neurites.time(i_nn) = time_frame;
					new_neurites.speed_el_tf(i_nn) = routes(best_idx).distance / (time_frame - previous_peaks.time(best_idx)); % Speed in electrodes per timeframe
					new_neurites.length_el(i_nn) = routes(best_idx).distance;
					new_neurites.speed_ms(i_nn) = new_neurites.speed_el_tf(i_nn) * (17.5/50); % pitch (17.5 um) / duration of time frame (50 us)
					new_neurites.amplitude(i_nn) = current_peaks.amplitude(j);
					new_neurites.voltage(i_nn) = current_peaks.voltage(j);
					new_neurites.peak_id_from(i_nn) = previous_peaks.id(best_idx);
					new_neurites.peak_id_to(i_nn) = current_peaks.id(j); % old: max(previous_peaks(:,'id')) + j;
					new_neurites.track_type(i_nn) = 2; % 2 = dynamic programming
					new_neurites.previous_peak_time(i_nn) = previous_peaks.time(best_idx);
				end
			end
			
			ret = new_neurites;
		end
		
		function ret = FilterAdaptive(obj, dt)
% 			ret = imgaussfilt3(dt, 0.7);
% 			
% 			return;
			obj.Output("FilterAdaptive()", 99);
			
			% Adaptive signal filtering (Wiener)
			time_points = 1:size(dt,3);
			dt_filtered = zeros(size(dt));
			for tp = time_points
				dt_filtered(:,:,tp) = wiener2(dt(:,:,tp),[obj.GetSettings().noise_filter.adaptive.size_1, obj.GetSettings().noise_filter.adaptive.size_2]);
			end
			ret = dt_filtered;
		end

		function ret = FilterImgaussfilt3(obj, dt)
			obj.Output("FilterImgaussfilt3()", 99);
			
			ret = imgaussfilt3(dt, obj.GetSettings().noise_filter.imgaussfilt3.sigma);
		end

		function ret = FilterImgaussfilt3_cort(obj, dt)
			obj.Output("FilterImgaussfilt3()", 99);
			
			ret = imgaussfilt3(dt, 0.7);
		end
		
		function ret = FilterBilateral(obj, dt)
			obj.Output("FilterBilateral()", 99);
			% Bilateral filtering
			% Requires a patch of the image which is used for smoothing (a patch of
			% noise in our case)
			time_points = 1:size(dt,3); % Number of measurements (time points in our case)
			
			filter_settings.intensity = 4;
			filter_settings.spatial_sigma = 7;

			dt_filtered = zeros(size(dt));
			for tp = time_points
				% Convert to L*a*b colorspace % Not sure we need this since our
				% data is not in RGB
				% rgb2lab() 

				% Select patch
		%         %test:
		%         figure();
		%         subplot(2,1,1)
		%         imshow(dt(150).signal)
		%         subplot(2,1,2)
		%         patch = imcrop(dt(150).signal, [180,80, 40, 40]); % Top right, the array is: [x, y, width, height]
		%         patch_std = std2(patch)^2;
		%         degree_smoothing = 2*patch_std;
		%         j = imbilatfilt(dt(150).signal, degree_smoothing, settings.spatial_sigma);
		%         imshow(j);
				% Select patch for smoothing
				patch = imcrop(dt(:,:,tp), [180,80, 40, 40]); % Top right, the array is: [x, y, width, height]
				patch_std = std2(patch)^2;
				degree_smoothing = 2 * patch_std * filter_settings.intensity;
				dt_filtered(:,:,tp) = imbilatfilt(dt(:,:,tp), degree_smoothing, filter_settings.spatial_sigma); % Filter using patch data, spatial sigma = bigger is "distant neighboring pixels contribute more to the Gaussian smoothing kernel"
			end

			%dt_filtered = dt_raw(3:end-2,3:end-2,:); % Make same size as the (kalman) filtered data
			ret = dt_filtered;
		end

		
		function nul = FindPeaks(obj)
			obj.Output("FindPeaks()", 0);
			
			obj.CalculateAbsoluteNoiseThresholds();
			
			time_frames = obj.GetSignalInit().time : size(obj.GetDtTracking(),3);
			time_frames_backwards = flip(time_frames);
			
			obj.AddPeaks(obj.GetSignalInit());
			
			direction = "backward";
			for time_frame = time_frames_backwards
				obj.AddPeaks(obj.FindPeaksInTimeFrame(time_frame, 'high',direction));
				obj.AddPeaks(obj.FindPeaksInTimeFrame(time_frame, 'medium',direction));
				obj.AddPeaks(obj.FindPeaksInTimeFrame(time_frame, 'low',direction));
				
				if obj.IsTooNoisy(); break; end
			end
			
			direction = "forward";
			for time_frame = time_frames
				obj.AddPeaks(obj.FindPeaksInTimeFrame(time_frame, 'high',direction));
				obj.AddPeaks(obj.FindPeaksInTimeFrame(time_frame, 'medium',direction));
				obj.AddPeaks(obj.FindPeaksInTimeFrame(time_frame, 'low',direction));
				
				if obj.IsTooNoisy(); break; end
			end
		end
		
		function ret = FindPeaksInTimeFrame(obj, time_frame, noise_threshold_level, direction)			
			ret = [];
			
			dt_voltage = obj.GetDtVoltage();
			dt_amplitude = obj.GetDtAmplitude();
			dt_tracking = obj.GetDtTracking();
			
			% Create local variables to improve readability
			dt_voltage_tf = dt_voltage(:,:,time_frame);
			dt_amplitude_tf = dt_amplitude(:,:,time_frame);
			dt_tracking_tf = dt_tracking(:,:,time_frame);
			signal_init = obj.GetSignalInit();
			settings = obj.GetSettings();
			new_peaks = obj.CreateEmptyTableCopy(obj.GetPeaks(), 0);
			
			% Set all values below the noise threshold to 0 (so we don't find peaks in that data)
			electrodes_x = size(obj.GetDtTracking(),2);
			electrodes_y = size(obj.GetDtTracking(),1);
			[X, Y] = meshgrid(1:electrodes_x,1:electrodes_y);
			
			noise_threshold = (settings.noise_threshold.(obj.GetNeuronType()).absolute.(noise_threshold_level));
			
			d_cd_f_thresholded = dt_tracking_tf;
			d_cd_f_thresholded(d_cd_f_thresholded < noise_threshold) = 0;
			max_tmp = imregionalmax(d_cd_f_thresholded, settings.peak_conn); % <= This finds peaks for us
    
			% If everything is the same value (e.g. 0 / no signal) everything is considered a peak by imregionalmax() (but actually there are no peaks)
			if sum(max_tmp(:)) == electrodes_x * electrodes_y
				return;
			end
			
			peaks_x = X(max_tmp);
			peaks_y = Y(max_tmp);

			% Activation radius: Peaks must be within a reasonable distance of the initiation site
			peaks_x = peaks_x(abs(X(max_tmp)-signal_init.x)<((time_frame-signal_init.time)*settings.max_speed) & abs(Y(max_tmp)-signal_init.y)<((time_frame-signal_init.time)*settings.max_speed));
			peaks_y = peaks_y(abs(X(max_tmp)-signal_init.x)<((time_frame-signal_init.time)*settings.max_speed) & abs(Y(max_tmp)-signal_init.y)<((time_frame-signal_init.time)*settings.max_speed));
			peaks_x_y = [peaks_x, peaks_y];
			 
			% Low/medium threshold peaks must be near previously detected peaks
			if settings.neighbor_dist.(noise_threshold_level) > 0
				
				peaks = obj.GetPeaks();
				% Get previous peaks, so we can see if the current peaks are close to them
				if direction == "forward"
					previous_peaks = peaks(peaks.time >= time_frame-settings.neighbor_lookback_time & peaks.time < time_frame,:);
				elseif direction == "backward"
					previous_peaks = peaks(peaks.time <= time_frame+settings.neighbor_lookback_time & peaks.time > time_frame,:);
				end
				
				% If there are not previous peaks, all peaks we found are invalid
				if isempty(previous_peaks)
					peaks_x_y = [];
				end
				
				% Check distance of current peak to previous peak, remove if too far away
				remove_indexes = [];
				for i_pxy = 1:size(peaks_x_y,1)					
					distances = sqrt(sum(    ( peaks_x_y(i_pxy,:) - previous_peaks{:,{'x','y'}} )     .^2, 2));
					if(min(distances) > settings.neighbor_dist.(noise_threshold_level))
						remove_indexes = [remove_indexes, i_pxy]; % Remove peak if it's too far away from its neighbors
					end
				end
				peaks_x_y(remove_indexes,:) = []; 
			end
    
			% Improve location accuracy and return
			if ~isempty(peaks_x_y)
				[peaks_x_improved, peaks_y_improved] = obj.ImproveAccuracyPeaks(dt_tracking_tf, peaks_x_y);
				
				% Make a local empty copy of the peaks table				
				new_peaks = obj.CreateEmptyTableCopy(obj.GetPeaks(), size(peaks_x_improved,1));
				new_peaks.x = peaks_x_improved(:,1);
				new_peaks.y = peaks_y_improved(:,1);
				new_peaks.time = zeros(size(peaks_x_improved,1),1)+time_frame;
				new_peaks.voltage = dt_voltage_tf(sub2ind(size(dt_voltage_tf), peaks_x_y(:,2), peaks_x_y(:,1)));
				new_peaks.amplitude = dt_amplitude_tf(sub2ind(size(dt_amplitude_tf), peaks_x_y(:,2), peaks_x_y(:,1)));
				new_peaks.confidence = zeros(size(peaks_x_improved,1),1); % TODO
				new_peaks.is_signal_init = zeros(size(peaks_x_improved,1),1);
			end
			ret = new_peaks;
		end
		
		function [peaks_x_improved, peaks_y_improved] = ImproveAccuracyPeaks(obj, dt_tracking_tf, peaks_x_y)
			peaks_x_improved = [];
			peaks_y_improved = [];

			if isempty(peaks_x_y)
				return;
			end

			electrodes_x = size(dt_tracking_tf,2);
			electrodes_y = size(dt_tracking_tf,1);

			peaks_x = peaks_x_y(:,1);
			peaks_y = peaks_x_y(:,2);
			peaks_x_improved = peaks_x;
			peaks_y_improved = peaks_y;
			for j=1:size(peaks_x,1)
				b = 0;
				if peaks_x(j,:) > 1 && peaks_x(j,:) < electrodes_x % Skip edges because we need the surrounding values
					a = dt_tracking_tf(peaks_y(j,:), peaks_x(j,:)-1:peaks_x(j,:)+1);
					b = (1-(max(a) - a(3)) / (max(a)-min(a)))*0.5;
					b = b + (-1-(a(1) - max(a)) / (max(a)-min(a)))*0.5;
				end
				peaks_x_improved(j,:) = peaks_x(j,:)+b;
			end
			for j=1:size(peaks_y,1)
				b = 0;
				if peaks_y(j,:) > 1 && peaks_y(j,:) < electrodes_y % Skip edges because we need the surrounding values
					a = (dt_tracking_tf(peaks_y(j,:)-1:peaks_y(j,:)+1, peaks_x(j,:)))';
					b = (1-(max(a) - a(3)) / (max(a)-min(a)))*0.5;
					b = b + (-1-(a(1) - max(a)) / (max(a)-min(a)))*0.5;
				end
				peaks_y_improved(j,:) = peaks_y(j,:)+b;
			end
		end
		
		function new = LeftMergeStruct(obj, orig, repl)
			% Merge repl into orig
			% Values in repl will overwrite values in orig
			% orig values that are not in repl will stay the same as orig
			% repl values that are not in orig will be added to the struct
			fieldnames_repl = fieldnames(repl);
			for i_fn = 1:size(fieldnames_repl, 1)
				if isstruct(repl.(fieldnames_repl{i_fn}))
					o = orig.(fieldnames_repl{i_fn});
					r = repl.(fieldnames_repl{i_fn});
					orig.(fieldnames_repl{i_fn}) = obj.LeftMergeStruct(o,r);
				else
					orig.(fieldnames_repl{i_fn}) = repl.(fieldnames_repl{i_fn});
				end
			end
			new = orig;
		end

		function nul = LoadExportedData(obj)
			obj.Output("LoadExportedData()", 0);
			
			file_base_name = obj.GetRecordingId().path_tracking;
			
			% If the table has more columns compared to the saved data (maybe we add new data over time)
			% then create a table with the newest format, and fill it with the saved data
			neurites = obj.CreateEmptyTableCopy(obj.GetNeurites(),0);
			neurites_loaded = load(strcat(file_base_name,'.neurites.mat'));
			neurites_loaded = neurites_loaded.neurites;
			neurites(1:size(neurites_loaded,1),neurites_loaded.Properties.VariableNames) = neurites_loaded(:,:);
			obj.neurites = neurites;

			% If the table has more columns compared to the saved data (maybe we add new data over time)
			% then create a table with the newest format, and fill it with the saved data
			peaks = obj.CreateEmptyTableCopy(obj.GetPeaks(),0);
			peaks_loaded = load(strcat(file_base_name,'.peaks.mat'));
			peaks_loaded = peaks_loaded.peaks;
% 			peaks_loaded
			if ~istable(peaks_loaded) % Using the old matrix format. Convert it to table.
				peaks_loaded(:,7) = zeros(1,size(peaks_loaded,1)); % voltage amplitudes
				% get amplitudes for the peaks
				dt_amplitude = obj.GetDtAmplitude();
				dt_voltage = obj.GetDtVoltage();
				for i_p = 1:size(peaks_loaded,1)
					x = peaks_loaded(i_p,1);
					y = peaks_loaded(i_p,2);
					t = peaks_loaded(i_p,3);
					peaks_loaded(i_p, 4) = dt_voltage(y,x,t);
					peaks_loaded(i_p, 7) = dt_amplitude(y,x,t);
				end
				peaks{1:size(peaks_loaded,1), {'x','y','time', 'voltage', 'id','confidence','amplitude'}} = peaks_loaded(:,:);
				
				peaks{peaks.amplitude==max(peaks.amplitude),{'is_signal_init'}} = 1;
			else
				peaks(1:size(peaks_loaded,1),peaks_loaded.Properties.VariableNames) = peaks_loaded(:,:);
			end
			obj.peaks = peaks;
			
			% Load settings. Default values will be used for values that are missing from the loaded data (might occur if we added new settings over time)
			obj.settings=[];
			default_settings = obj.GetSettings();
			loaded_settings = load(strcat(file_base_name,'.settings.mat'));
			loaded_settings = loaded_settings.settings;
			obj.settings = obj.LeftMergeStruct(default_settings, loaded_settings);
			
			dt_tracking = load(strcat(file_base_name,'.dt_filtered.mat'));
			obj.dt_tracking = dt_tracking.dt_filtered;
		end
		
		function nul = Output(obj, text, verbose_level)
			if verbose_level >= 1
				disp(text)
			end
		end
		
		function nul = RemovePrePostPeaks(obj)
			% Due to the shape of the signal, a small peak before or after the main peak can be detected
			% We want to remove these, otherwise we will get overlapping "double" tracks
			peaks = obj.GetPeaks();
			max_iter = 99999;
			cur_iter = 1;
			peak_ids_to_remove = nan(max_iter,1);	% Pre-allocate the memory, otherwise we'll get memory leaks. See: https://se.mathworks.com/matlabcentral/answers/21004-is-there-an-elegant-way-to-create-dynamic-array-in-matlab
			for i_p = 1:size(peaks,1)
				
				% if the time between the peaks is between 4 and 10, this is likely a pre or post peak
				peaks_close_time_range = peaks(abs(peaks.time - peaks.time(i_p)) > 6 & abs(peaks.time - peaks.time(i_p)) < 13,:);
				
				if isempty(peaks_close_time_range); continue; end
				
				% Get closest peak
				distances =  sqrt(sum((peaks_close_time_range{:,{'x','y'}} - peaks{i_p,{'x','y'}}).^2,2));
				
				if isempty(distances(distances < 3)); continue; end
				
				close_peaks = peaks_close_time_range(distances < 3,:);
				
				% Pre and post peaks should have a lower amplitude than 70% the main peak.
				if peaks.amplitude(i_p) < (max(close_peaks.amplitude) * 0.70)
					% Delete
					peak_ids_to_remove(cur_iter) = peaks.id(i_p);
					cur_iter = cur_iter + 1;
					if cur_iter > max_iter
						obj.Output("Reached max nr of peaks to remove.", 99);
						break;
					end
				end
			end
			peak_ids_to_remove = peak_ids_to_remove(~isnan(peak_ids_to_remove),:);
			obj.DelPeaks(peak_ids_to_remove);
		end
		
		
		% CHECKS / VALIDATION
		function valid = IsValidRecordingId(obj)
			if isempty(obj.GetRecordingId())
				disp("ERROR: Invalid recording id provided. Set recording ID using NewTrack() or LoadTrack()");
				valid = 0;
				return
			end
			if size(obj.GetRecordingId(),1) > 1
				disp("ERROR: Provide only one recording ID");
				valid = 0;
				return
			end
			valid = 1;
		end
		
		function valid = IsToBeAxtrackted(obj)
			valid = 1;
			if ismember(obj.recording_id.recording_id , obj.GetAxdb().recording_id)
				valid = 0;
			end
		end
		
		function valid = IsTooNoisy(obj)
			% TODO
			valid = 0;
		end
		
		function valid = HasValidCelltype(obj)
			valid = 1;
			if obj.GetNeuronType() == obj.neuron_types(3) % unknown cell type
				valid = 0;
			end
		end
	end
end

































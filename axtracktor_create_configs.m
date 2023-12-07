%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Divide neuron into N electrode batches so we can track individual
% branches.
%  neurites: table containing neurite information
%  dirname: output directory name
%  dfile: output file name
%  settings: various settings. e.g. initation start time and such
% Return: Nothing
% 
% Author: Markus de Ruijter; Punga Lab, Uppsala University; 2019

function null = axtracktor_create_configs(neurites2, dirname, dfile, settings)

    disp('Clustering neurites..');

    % So now we just identify segments that correspond to X/3 electrodes
    length_per_cluster = 200; % * 3 = 600 electrodes

    % Iterate from the highest branch level and move down until all points have
    % a cluster assigned
    % First make small clusters, individual branches. Then merge/optimize into
    % larger clusters
    cluster_idx = 1;
    % clusters = dataset([],[],[],[],'VarNames',{'neurite_idx','total_length','branching_segment','merged'}); % unused but this is how you create and empty dataset
    clusters = struct();
    neurites2 = sortrows(neurites2, {'branch_level','segment_id'});
    neurite_idx = find(neurites2.branch_level == 0)';
    clustered_neurite_idx = [];
    we_need_to_find_more_clusters = 1;
    max_branch_level = max(neurites2.branch_level);
    while we_need_to_find_more_clusters 
        new_neurite_idx = [];
        for i = 1:size(neurites2,1)
            neurite = neurites2(i,:);
            if ismember(i, clustered_neurite_idx) 
                continue;
            end
            if neurite.branch_length < length_per_cluster
                %neurites2.cluster_idx(find(neurites2.neurite_id >= neurite.neurite_id & neurites2.neurite_id < (neurite.neurite_id + 1 * 10^(max_branch_level - neurite.branch_level)))') = cluster_idx;


                % Get all indexes in this branch, including the path to the
                % initiation site
                branch_neurite_idx = find(neurites2.neurite_id >= neurite.neurite_id... 
                    & neurites2.neurite_id < (neurite.neurite_id + 1 * 10^(max_branch_level - neurite.branch_level))...
                    & neurites2.line_id >= -1)'; % Line id -1 = pseudo connection

                for bl = 1:neurite.branch_level
    %                 disp(bl)
    %                 disp(neurite.neurite_id)
    %                 disp(roundn(neurite.neurite_id, max_branch_level-(bl-1)))
                    branch_neurite_idx = [branch_neurite_idx, find(neurites2.neurite_id == extractBefore(neurite.neurite_id, bl+1)...
                        & neurites2.track_type ~= 3)']; % track_type 3 = pseudo connection
                end
                clusters(cluster_idx).neurite_idx = branch_neurite_idx;
                clustered_neurite_idx = [clustered_neurite_idx, branch_neurite_idx];
                clusters(cluster_idx).total_length = sum(neurites2.length(branch_neurite_idx));
                clusters(cluster_idx).branching_segment = i;
                clusters(cluster_idx).merged = 0;
                cluster_idx = cluster_idx + 1;
            end

        end

        neurite_idx = new_neurite_idx;

        if isempty(neurite_idx)
            we_need_to_find_more_clusters = 0;
        end
    end

    % Merge clusters
    for i = 1:size(clusters,2)
        best_merge_idx = 0;
        best_merge_diff = -1;
        for j = (i+1):size(clusters,2)
            cluster_neurite_idx = unique([clusters(i).neurite_idx, clusters(j).neurite_idx]);
            size_both = clusters(i).total_length + clusters(j).total_length;
            size_merged = sum(neurites2.length(cluster_neurite_idx));
            size_diff = size_both - size_merged;
            if size_diff > best_merge_diff && size_merged < length_per_cluster
                best_merge_idx = j;
                best_merge_diff = size_diff;
            end
        end

        if best_merge_idx > 0
            clusters(best_merge_idx).neurite_idx = unique([clusters(best_merge_idx).neurite_idx, clusters(i).neurite_idx]);
            clusters(best_merge_idx).total_length = sum(neurites2.length(clusters(best_merge_idx).neurite_idx));
            clusters(i).merged = 1;
        end
    end

    clusters = clusters(~[clusters.merged]);

    for i = 1:size(clusters,2)
        cluster_electrodes = round(unique([neurites2.x_from(clusters(i).neurite_idx), neurites2.y_to(clusters(i).neurite_idx)], 'rows'));


        configuration_area = [[cluster_electrodes(:,1)-1; ... % X coords
            cluster_electrodes(:,1)-1;...
            cluster_electrodes(:,1)-1;...
            cluster_electrodes(:,1);...
            cluster_electrodes(:,1);...
            cluster_electrodes(:,1);...
            cluster_electrodes(:,1)+1;...
            cluster_electrodes(:,1)+1;...
            cluster_electrodes(:,1)+1],...
            [cluster_electrodes(:,2)-1;... % Y coords
            cluster_electrodes(:,2);...
            cluster_electrodes(:,2)+1;...
            cluster_electrodes(:,2)-1;...
            cluster_electrodes(:,2);...
            cluster_electrodes(:,2)+1;...
            cluster_electrodes(:,2)-1;...
            cluster_electrodes(:,2);...
            cluster_electrodes(:,2)+1]];
        configuration_area = unique(configuration_area, 'rows');
        clusters(i).configuration_area = configuration_area;
        clusters(i).config_area_size = size(configuration_area,1);

    end
    
    % Save the data
    save_path = './output/data/';
    mkdir(save_path);
    save(strcat(save_path,dirname,'_',dfile, settings.file_output_postfix,'.configs.mat'), 'clusters');

    % Use subtightplot 
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.025]);
    if ~make_it_tight,  clear subplot;  end
    
    f = figure('color','black', 'Visible','off');
    subplot(1,1,1);
    whitebg('black');   % Change color of area behind/around the figure
    scatter(settings.signal_init.x, settings.signal_init.y, 500, 'white', 'filled','o');%,cls,'filled','s');  
    alpha(0.1);
    hold on;

    color_palette = hsv(max([size(clusters,2), 4]));
    total_config = [];

    clearvars fake_plot;
    for i = 1:size(clusters,2)
        cluster_config = clusters(i).configuration_area;

        cluster_config(:,[3,4,5]) = repelem(color_palette(i,:), size(cluster_config,1), 1);

        total_config = [total_config; cluster_config];
        fake_plot(i) = scatter(nan,nan,16, color_palette(i,:),'filled','o');
        hold on;
    end

    total_config_ds = dataset({total_config,'x','y','r','g','b'});
    total_config_ds_grp = grpstats(total_config_ds,{'x','y'});

    scatter(total_config_ds_grp.x, total_config_ds_grp.y, total_config_ds_grp.GroupCount*8, [total_config_ds_grp.mean_r, total_config_ds_grp.mean_g, total_config_ds_grp.mean_b], 'filled','o');
    alpha(0.8);
    hold on;
    lines_plot = line([neurites2.x_from'; neurites2.x_to'], [neurites2.y_from'; neurites2.y_to']);
    hold on; 
    daspect([1 1 1]);

    color_palette = jet(size(unique((neurites2.line_id)),1));
    [null, clr_ids] = ismember(neurites2.line_id', unique((neurites2.line_id)));
    line_cls2 = color_palette(clr_ids,:);
    for k = 1 : length(lines_plot)
        %set(lines_plot(k), {'Color','LineWidth'}, {line_cls2(k,:), 3});
        set(lines_plot(k), {'Color','LineWidth'}, {'white', 2});
    end
    leg = legend(fake_plot, num2str([clusters.config_area_size]'));
    title(leg,'# electr.');
    
    f.InvertHardcopy = 'off'; % Preserve colors when saving
    rez=200; %resolution (dpi) of final graphic
    figpos=getpixelposition(f); %dont need to change anything here
    resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
    set(f,'paperunits','points','papersize',figpos(3:4)/resolution,...
    'paperposition',[0 0 600 600]); %paperposition can be used to actually increase figure size

    save_path = './output/config_figs/';
    mkdir(save_path);
    print(f,fullfile(strcat(save_path,dirname,'_',dfile, settings.file_output_postfix,'.png')),'-dpng',['-r',num2str(rez)],'-opengl'); %save file
end































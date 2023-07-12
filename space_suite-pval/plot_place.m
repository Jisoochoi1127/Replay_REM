function plot_place(place_field_data, ms, interpol_behav, behav);
    
    place_field=place_field_data.PlaceField;
    interpolated_Y=interpol_behav.interpolated_Y;
    interpolated_X=interpol_behav.interpolated_X;
    interpolated_speed=interpol_behav.interpolated_speed;

    %% Plotting the results
    place_field_mask_img = place_field/max(place_field(:))*150;
    place_field_mask_img(:,:,2) = place_field/max(place_field(:))*150;
    place_field_mask_img(:,:,3) = 0;

    place_field_mask_img = imresize(place_field_mask_img,size(behav.background,1)./size(place_field_mask_img,1),'bilinear');
    place_field_mask_img = uint8(place_field_mask_img);

    color_vector(:,1) = linspace(0,1,length(ms.time));
    color_vector(:,3) = color_vector(:,1);
    color_vector(:,2) = flipud(color_vector(:,1));

    % Converting back to pixels
    interpolated_X_pix = interpolated_X*behav.ROI(3)/behav.trackLength;
    interpolated_Y_pix = interpolated_Y*behav.ROI(3)/behav.trackLength;

    RUN_binarized_trace = ms.Binary(:,place_field_data.Cell_ID);
    RUN_binarized_trace(interpolated_speed<place_field_data.Params.min_speed) = 0; % Remove non-running activity

    %% Spatial firing
    clf
    plotting_fig = gcf;
    set(plotting_fig,'Name',strcat('Cell ID: ',num2str(place_field_data.Cell_ID)),'NumberTitle','off')
    subplot(7,4,[9 10 13 14]);
    if isfield(behav,'background')
    imshow(3*behav.background); hold on
    end
    plot(interpolated_X_pix,interpolated_Y_pix,'color',[1 1 1 0.3]); hold on;
    scatter(interpolated_X_pix(logical(RUN_binarized_trace)),interpolated_Y_pix(logical(RUN_binarized_trace)), [100],color_vector(logical(RUN_binarized_trace),:), '.');
    pf_mask = imshow(place_field_mask_img);
    set(pf_mask, 'AlphaData', 0.5);
    ax1=gca;
    colormap(ax1,color_vector);
    %cb=colorbar(ax1,'southoutside','Ticks',[0,1],'TickLabels',{'Start','End'});
    %cb.Label.String = 'Time';
    daspect([1 1 1])
    ax1.YDir = 'Reverse';
    set(ax1, 'Visible','off');hold off
    title 'Activity across time/space'
    cb=colorbar(ax1,'EastOutside','Ticks',[0,1],'TickLabels',{'Start','End'});

    % Total activity      
     subplot(7,4,[1 2 5 6]);
    pz=zeros(3,length(place_field_data.spatial_matrix));
    pz(2,:)=place_field_data.spatial_matrix;
    pz=imgaussfilt(pz);
    pcolor(pz);
    shading interp
    hold on
    title('Unmasked firing map')
    hcb=colorbar;
    title(hcb,'Calcium activity')

    shading interp;
    set (gca,'DataAspectRatio',[1 1 1],'YDir','Reverse',...
        'XTick',[],...
        'YTick',[]);
    
          
    subplot(7,4,[11 12 15 16]);
    pz=zeros(3,length(place_field_data.spatial_matrix));
    pz(2,:)=place_field_data.PlaceField;
    pz=imgaussfilt(pz);
    pcolor(pz);
    shading interp
    hold on
    title('Masked firing map')
    hcb=colorbar;
    title(hcb,'Calcium activity')
    if place_field_data.numPlaceFields > 0        
        for i=1:length(place_field_data.PlaceFieldCentroid)  
            scatter(place_field_data.PlaceFieldCentroid{1,i}(1),place_field_data.PlaceFieldCentroid{1,i}(2)*2, 200, [1,1,1], '+')
        end
    end

    shading interp;
    set (gca,'DataAspectRatio',[1 1 1],'YDir','Reverse',...
        'XTick',[],...
        'YTick',[]);

    subplot(7,4,[3 4 7 8]);
    imagesc(place_field_data.lap_mat);
    ylabel 'Laps'
    xlabel 'Bin #'    
    title 'Trajectories'    
    
    hcb=colorbar;
    title(hcb,'Calcium activity')
    


    subplot(7,4,[17:20]);
    plot(ms.time/1000,zscore(ms.RawTraces(:,place_field_data.Cell_ID)),'color','black');
    hold on
    plot(ms.time/1000,-ms.Binary(:,place_field_data.Cell_ID)+max(zscore(ms.RawTraces(:,place_field_data.Cell_ID)))+1,'color','red');
    ax7=gca;
    ax7.XLim= [0 ms.time(end)/1000];
    ax7.YLim= [min(zscore(ms.RawTraces(:,place_field_data.Cell_ID))) max(zscore(ms.RawTraces(:,place_field_data.Cell_ID))+0.5)];
    
    subplot(7,4,[21:24]);
    plot(behav.time/1000,behav.position);
    ax8=gca;
    ax8.XLim= [0 ms.time(end)/1000];
    ax8.YLim= [0 100];
    title 'Position (cm)'

    subplot(7,4,[25:28]);
    plot(behav.time/1000,behav.speed);
    ax9=gca;
    ax9.XLim= [0 ms.time(end)/1000];
    title 'Speed (cm.s-1)'    
    drawnow
end

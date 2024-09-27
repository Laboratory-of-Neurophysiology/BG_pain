    close all; clear all;
    
    time_frequency =  2 %Hz
    
    %High-pass filter (HF = 1/0)
    HF = 0
    % sheet
    sheet_NO =  1 
    file_name = '.xlsx' % write the file name
    
    low_value = 0.1 %(0-1.0: 0.3 default)
    threshold = 0.01 %normalized f_signal

    window_size = 100
        
    gaussian_window_size = 10 %(0-10)
    
    low_high_para = 'low' % low or high 
    Wn_para = 1 %cut-off value(1-5)
       
    %% Interval analysis
     N_compartment = 1;
        compart_info = zeros(N_compartment,2); 
               	compart_info(1,:) = [0 60];

    %% LOAD
      
    for i = 1:sheet_NO
        table(:,:,i) = xlsread(file_name, i);
    end

    table(:,1,:) = []; %time column

    [r, c, z] = size(table);
    sec_length = r/time_frequency;
    xval = 1/time_frequency:1/time_frequency:sec_length; 
    
    %% normalize
    f0_vector = zeros(z,c);
    f0 = zeros(r,c,z);   
    
    f_signal = zeros(r,c,z);
    
    for i = 1:z
        table_sorted = sort(table(:,:,i));
        f0_vector(i,:) = mean(table_sorted(1:round(r*low_value),:));
        for j = 1:r
            f_signal(j,:,i) = (table(j,:,i) - f0_vector(i,:))./f0_vector(i,:);
        end
        
        f_signal_sorted = sort(f_signal(:,:,i));
        f0_signal_std(i,:) = std(f_signal_sorted(1:round(r*low_value),:));
        f0_std(:,:,i) = repmat(f0_signal_std(i,:),r,1);
    end
       
    for i = 1:z
        figure('position', [100, 100, 1000, 200])
        plot(xval, f_signal(:,:,i))
              xlim([0 30])
        title(i); xlabel('Time (sec)');
    end
        
    %f_signal display
    f_display = zeros(r,c,z);
    for i = 1:z
        for j = 1:c
            f_display(:,j,i) = f_signal(:,j,i) -m*(j-1);
        end
    end

    for i = 1:z
        figure('position', [100, 100, 1000, 200])
        plot(xval, f_display(:,:,i))
%         ylim([-0.05 0.2])
        title(i); xlabel('Time (sec)');
    end
        
    %% frequency filtering
    %High-pass/low-pass filter
    %filter design
    if HF == 1;
       
        Fs = 100; %sampling rate 
        n = 10; %order 5-10 recommended
        Wn = Wn_para; %cut-off frequency (½Â¾ð PC->2, BG->?)
        Fn = Fs/2;
        [b a] = butter(n,Wn/Fn, low_high_para); %movement correction/PC-> high, BG->low 
    
        %filter 
        for i = 1:sheet_NO
            f_signal(:,:,i) = filtfilt(b,a,f_signal(:,:,i));
        end
    end
    
    for i = 1:z
        figure        
        plot(xval, f_signal(:,:,i))
        title(i); xlabel('Time (sec)');
    end
        
    %% gaussian filtering
   
    gaussFilter = gausswin(gaussian_window_size);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
    
    for i = 1:z
        shift_index = zeros(c,1);
        for j = 1:c
            
            valid_trace = conv(f_signal(:,j,i),gaussFilter, 'valid');
            full_trace = conv(f_signal(:,j,i),gaussFilter, 'full');
            
            f_signal_z_gaussian(:,j,i) = valid_trace;
                       
            cir_matrix = zeros(length(full_trace), length(full_trace));
            for k=1:length(full_trace)
                cir_matrix(:,k) = circshift(full_trace, k-1);
            end
            
            ttt = corr(cir_matrix(round(r/8):r-round(r/8),:), valid_trace(round(r/8):r-round(r/8)));
            temp = find(ttt ==  max(ttt));
            shift_index(j) = length(full_trace) - temp + 1;
            shift_index(j) % before xval
        end
            remove_fore = median(shift_index); % before xval
            remove_back = length(full_trace) - length(valid_trace) - remove_fore;
            
            
   end
   new_xval = xval(remove_fore+1: end); % new x axis after filtering
   
   %% update after new_xval
       for i =1:N_compartment
        for j =1:2
            [cc index] = min(abs(new_xval - compart_info(i,j)));
            compart_info(i,j) =  new_xval(index);
        end
       end
    %%
    
    %figure; plot(table_filtered)
    f_signal = f_signal_z_gaussian;
    [r c z] = size(f_signal);
    
    for i = 1:z
        figure     
        plot(new_xval, f_signal(:,:,i))
        title(i); xlabel('Time (sec)');
        
        figure('position', [100, 100, 1200, 300]);
        imagesc(f_signal'*80, [0 150])
        xlim([1 300])
        colormap jet
        colorbar; 
    end
    
    %% correlationi matrix 
    % + average correlation 
    average_correlation = zeros(z,1);

    %check the minimum correlation values across the sheets (min_corr)
    min_corr = 1;
    max_corr = 0;
    for i = 1:z
        min_corr = min([min(min(corrcoef(f_signal(:,:,i)))), min_corr]);

        ascend_sort  = sort(corrcoef(f_signal(:,:,i)), 'descend');
        max_corr = max(max(ascend_sort(2,:)), max_corr);
     end

    figure
    for i = 1:z
        % cell by cell correlation
        cell_correlation = corrcoef(f_signal(:,:,i));

        subplot(2,z,i);
        imagesc(cell_correlation,[min_corr max_corr])
        title(i)
        average_correlation(i) = mean(mean(cell_correlation));

        % time by time correlation
        time_correlation = corrcoef(f_signal(:,:,i)');

        subplot(2,z,i+z);
        imagesc(time_correlation)
        title(i)
    end

    figure
    plot(average_correlation,'ro--');
    title('average correlation')
    
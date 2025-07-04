function filled_final = detect_etalon_band(wt_671, f_671, wt_785, f_785, wl, f_axis, loop_times, visible)
% DETECT_ETALON_BAND Detect etalon fringe band from wavelet spectra
%
%   filled_final = detect_etalon_band(wt_671, wt_785, wl, f_axis, loop_times)
%
%   Inputs:
%       wt_671     - Wavelet coefficients of 671 nm ICA component (freq x wl)
%       f_671      - Frequency vector corresponding to wavelet coefficients of 671 nm (Hz)
%       wt_785     - Wavelet coefficients of 785 nm ICA component (freq x wl)
%       f_785      - Frequency vector corresponding to wavelet coefficients of 785 nm (Hz)
%       wl         - Wavelength vector (1 x N)
%       f_axis     - Frequency vector from CWT (1 x F)
%       loop_times - Integer indicating iteration (affects threshold)
%       visible    - whether to display intermediate processing figures (logical flag)
%
%   Output:
%       filled_final - Binary mask (freq x wl) indicating etalon band region
    
    if nargin < 8
        visible = false;  % 设置默认值
    end

    % Normalize wavelet coefficients
    gray_671 = abs(wt_671) / max(abs(wt_671(:)));
    gray_785 = abs(wt_785) / max(abs(wt_785(:)));

    % Canny edge detection with adaptive threshold using loop control
    switch loop_times
        case 1
            can_671 = 0.04;
            can_785 = 0.04;
        case 2
            can_671 = 0.03;
            can_785 = 0.03;
        otherwise
            can_671 = 0.005;
            can_785 = 0.005;
    end

    edges_671 = edge(gray_671, 'Canny', can_671);
    edges_785 = edge(gray_785, 'Canny', can_785);

    %% Visulization
    if visible
        figure;
        subplot(2,2,2);
        % imshow(edges_671);
        surface(wl,f_671,edges_671);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colorbar;
        colormap(gca,'Gray');
        title('Edge of Frequency-Wavelength spectrum (671nm)')
    
        subplot(2,2,1);
        surface(wl,f_671,gray_671);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colorbar;
        colormap(gca,'jet');
        caxis([0 0.05]);
        title('Frequency-Wavelength spectrum of etalon ic (671 nm)')
    
        subplot(2,2,4);
        % imshow(edges_671);
        surface(wl,f_785,edges_785);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colorbar;
        colormap(gca,'Gray');
        title('Edge of Frequency-Wavelength spectrum (785nm)')
    
        subplot(2,2,3);
        surface(wl,f_785,gray_785);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colorbar;
        colormap(gca,'jet');
        caxis([0 0.05]);
        title('Frequency-Wavelength spectrum of etalon ic (785 nm)')
    end
%%

    % Cross-correlation of edge maps
    C = double(edges_671) .* double(edges_785);

    % Region limits (based on wavelength & frequency)
    switch loop_times
        case 1
            w_min = 860;
            w_max = 960;
            f_min = 0.2;
            f_max = 0.8;
        case 2
            w_min = 890;
            w_max = 960;
            f_min = 0.25;
            f_max = 0.8;
        case 3
            w_min = 880;
            w_max = 960;
            f_min = 0.2;
            f_max = 0.8;
        case 4
            w_min = 900;
            w_max = 960;
            f_min = 0.2;
            f_max = 1.0;
        case 5
            w_min = 880;
            w_max = 910;
            f_min = 0.33;
            f_max = 1.0;
        case 6
            w_min = 886;
            w_max = 897;
            f_min = 0.5;
            f_max = 1.0;
        case 7
            w_min = 902;
            w_max = 960;
            f_min = 0.15;
            f_max = 0.45;
        otherwise
            w_min = 900;
            w_max = 960;
            f_min = 0.8;
            f_max = 0.8;
    end

    [~, w_min_pos] = min(abs(wl - w_min));
    [~, w_max_pos] = min(abs(wl - w_max));
    [~, f_min_pos] = min(abs(f_axis - f_min));
    [~, f_max_pos] = min(abs(f_axis - f_max));

    % Create binary mask
    mask = zeros(size(C));
    mask(f_max_pos:f_min_pos, w_min_pos:w_max_pos) = 1;
    C_cleaned = C .* mask;
    %%
    if visible
        figure;
        subplot(1,2,1);
        surface(wl,f_785,C);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colormap(gray);
        colorbar;
        xlabel('Wavelength (nm)');
        ylabel('Frequency');
        title('XCorr of noise content t-f');
    
        subplot(1,2,2);
        surface(wl,f_785,C_cleaned);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colormap(gray);
        colorbar;
        xlabel('Wavelength (nm)');
        ylabel('Frequency');
        title('XCorr of noise content t-f')
    end
    
    %% Extract top and bottom boundary
    [height, width] = size(C_cleaned);
    top_boundary = nan(1, width);
    bottom_boundary = nan(1, width);

    for col = 1:width
        col_data = C_cleaned(:, col);
        if sum(col_data) > 1
            top_idx = find(col_data, 1, 'first');
            bottom_idx = find(col_data, 1, 'last');
            if ~isempty(top_idx)
                top_boundary(col) = top_idx;
            end
            if ~isempty(bottom_idx)
                bottom_boundary(col) = bottom_idx;
            end
        end
    end

    % Detect invalid situation before interpolation
    if all(isnan(top_boundary)) || all(isnan(bottom_boundary))
        warning('No valid etalon edges found in loop %d, skipping.', loop_times);
        filled_final = [];  % 返回空矩阵，供主函数判断是否退出循环
        return;
    end

    x = 1:width;
    top_boundary_filled = top_boundary;
    bottom_boundary_filled = bottom_boundary;

    top_nan_idx = isnan(top_boundary); % nan出现的位置
    bottom_nan_idx = isnan(bottom_boundary);
    top_boundary_filled(top_nan_idx) = interp1(x(~top_nan_idx), top_boundary(~top_nan_idx), x(top_nan_idx), 'linear');
    bottom_boundary_filled(bottom_nan_idx) = interp1(x(~bottom_nan_idx), bottom_boundary(~bottom_nan_idx), x(bottom_nan_idx), 'linear');

    filled_img_top = zeros(height, width);
    filled_img_bottom = zeros(height, width);

    for col = 1:width
        if ~isnan(top_boundary_filled(col)) && ~isnan(bottom_boundary_filled(col)) && top_boundary_filled(col) ~= bottom_boundary_filled(col)
            top_idx = round(top_boundary_filled(col));
            bottom_idx = round(bottom_boundary_filled(col));
            filled_img_top(top_idx, col) = 1;
            filled_img_bottom(bottom_idx, col) = 1;
        elseif isnan(top_boundary_filled(col)) && isnan(bottom_boundary_filled(col))
            filled_img_top(:, col) = 0;
        elseif ~isnan(top_boundary_filled(col)) && ~isnan(bottom_boundary_filled(col)) && top_boundary_filled(col) == bottom_boundary_filled(col)
            filled_img_top(:, col) = 0;
        end
    end

    maxx = zeros(1, height);
    maxx_1 = zeros(1, height);

    for row = 1:height
        idx_top = find(filled_img_top(row, :));
        idx_bot = find(filled_img_bottom(row, :));

        if numel(idx_top) > 1
            maxx(row) = max(idx_top);
            first_nonzero_idx = find(maxx, 1, 'first');
            filled_img_top(row, min(idx_top):max(idx_top)) = 1;
            if row > 1 && min(idx_top) > max(maxx(1:row-1)) && first_nonzero_idx ~= row
                filled_img_top(row, max(maxx(1:row-1)):min(idx_top)) = 1;
            end
        end

        if numel(idx_bot) > 1
            maxx_1(row) = max(idx_bot);
            first_nonzero_idx = find(maxx_1, 1, 'first');
            filled_img_bottom(row, min(idx_bot):max(idx_bot)) = 1;
            if row > 1 && min(idx_bot) > max(maxx_1(1:row-1)) && first_nonzero_idx ~= row
                filled_img_bottom(row, max(maxx_1(1:row-1)):min(idx_bot)) = 1;
            end
        end
    end

    filled_final = filled_img_top + filled_img_bottom;
    for col = 1:width
        idx = find(filled_final(:, col));
        if numel(idx) > 1
            filled_final(min(idx):max(idx), col) = 1;
        end
    end

    %%
    if visible
        figure;
        subplot(2, 2, 1);
        surface(wl,f_785,C_cleaned);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colormap(gray);
        colorbar;
        xlabel('Wavelength (nm)');
        ylabel('Frequency');
        title('Original Binary Image');
    
        subplot(2, 2, 2);
        surface(wl,f_785,filled_img_top);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colormap(gray);
        colorbar;
        xlabel('Wavelength (nm)');
        ylabel('Frequency');
        title('Top Boundary Line');
    
        subplot(2, 2, 3);
        surface(wl,f_785,filled_img_bottom);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colormap(gray);
        colorbar;
        xlabel('Wavelength (nm)');
        ylabel('Frequency');
        title('Bottom Boundary Line');
    
        subplot(2, 2, 4);
        surface(wl,f_785,filled_final);
        shading interp;
        ylim([0 1]);
        xlim([810 960]);
        colormap(gray);
        colorbar;
        caxis([0 1]);
        xlabel('Wavelength (nm)');
        ylabel('Frequency');
        title('Wavelet Coefficient Mask');%% Reconstruct the noise by selecting valid frequency band
    end
%%

end

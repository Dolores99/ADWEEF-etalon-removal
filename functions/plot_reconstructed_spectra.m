function plot_reconstructed_spectra(recon_671, recon_785, wl_671, wl_785, depth_start)
% PLOT_RECONSTRUCTED_SPECTRA Visualize denoised Raman spectra (671 & 785 nm)
%
% Inputs:
%   recon_671    - M x N matrix of reconstructed spectra (671 nm)
%   recon_785    - M x N matrix of reconstructed spectra (785 nm)
%   wl_671       - 1 x N vector of wavelength values for 671 nm
%   wl_785       - 1 x N vector of wavelength values for 785 nm
%   depth_start  - starting row index for plotting (e.g., 2)

    figure2 = figure('Position', [2, 470, 638, 327]);
    colors = {'#0072BD','#EDB120','#77AC30','#A2142F','#7E2F8E','#D95319'};

    %% Plot 785 nm spectra
    ax1 = axes('Parent',figure2,'Position',[0.13 0.11 0.368125 0.815]);
    hold(ax1, 'on');
    for d = depth_start:(depth_start + 5)
        plot(wl_785, recon_785(d,:) + (d - depth_start) * 3000, 'Color', colors{d - depth_start + 1}, 'LineWidth', 2);
    end
    xlim([600 1800]);
    set(ax1, 'XTick', 600:200:2000, 'YTick', []);
    yticklabels([]);
    % 去掉右边框，但保留上侧框线
    set(axes1, 'Box', 'off');  % 关闭默认的边框
    set(axes1, 'YColor', 'k', 'YTick', []); % 保留左侧框线
    set(axes1, 'XColor', 'k'); % 保留下侧框线
    ylabel('Light Intensity (a.u.)','Color', 'K', 'FontName', 'Arial','FontWeight', 'bold');
    % 获取当前 Y 轴的范围
    yl = ylim;
    
    % 设置 Y 轴的最小值，并保持最大值自动调整
    ylim([-480 yl(2)]);  % 将最小值设置为 0，最大值保持原来的自动设置
    
    x_limits = axes1.XLim;
    y_limits = axes1.YLim;
    line([x_limits(1) x_limits(2)], [y_limits(2) y_limits(2)], 'Color', 'k');
    
    
    hold off;

    %% Plot 671 nm spectra
    ax2 = axes('Parent',figure2,...
    'Position',[0.520340909090909 0.11 0.29 0.815]);
    hold(ax2, 'on');
    for d = depth_start:(depth_start + 5)
        plot(wl_671, recon_671(d,:) + (d - depth_start) * 3000, 'Color', colors{d - depth_start + 1}, 'LineWidth', 2);
    end
    xlim([2800 3800]);
    set(ax2, 'XTick', [2800 3000 3200 3400 3600 3800]);
    yl = ylim;

    % 设置 Y 轴的最小值，并保持最大值自动调整
    ylim([-480 yl(2)]);  % 将最小值设置为 0，最大值保持原来的自动设置
    
    % 去掉左边框，但保留上侧框线
    set(axes2, 'Box', 'off');  % 关闭默认的边框
    set(axes2, 'YColor', 'none');  % 移除左侧框线
    set(axes2, 'XColor', 'k');  % 保留下侧框线
    
    x_limits = axes2.XLim;
    y_limits = axes2.YLim;
    line([x_limits(1) x_limits(2)], [y_limits(2) y_limits(2)], 'Color', 'k');
    line([x_limits(2) x_limits(2)], [y_limits(1) y_limits(2)], 'Color', 'k', 'LineWidth', 0.5);
    
    
    
    hold off;

end

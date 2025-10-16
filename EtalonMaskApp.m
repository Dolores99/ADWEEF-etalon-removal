function EtalonMaskApp(lambda, spectrum)
    % EtalonMaskApp 交互式在 CWT F–W 图上绘制掩膜，支持自由画笔、矩形和多边形，软抑制 etalon 分量并回看光谱
    % 使用: EtalonMaskApp(lambda, spectrum)
    % lambda: 1×N 波长 (nm)
    % spectrum: 1×N 对应强度（原始或 etalon-like IC）
    
    % 输入参数检查，确保是向量，并转为行向量 (方便后续矩阵运算)
    assert(isvector(lambda) && isvector(spectrum), 'lambda & spectrum must be vectors.');
    lambda = lambda(:).'; spectrum = spectrum(:).';
    N = numel(lambda); % 数据点数量
    
    % ---------------- UI (用户界面) ----------------
    % 创建主 App 窗口
    ui = uifigure('Name','Etalon Masking (F–W Interactive)','Position',[100 100 1000 720]);
    
    % 创建网格布局管理器
    g = uigridlayout(ui,[6 9]);
    % 定义行高：30(参数) 30(参数) 1x(轴区) 30(按钮) 1x(轴区) 30(指标)
    g.RowHeight = {50 30 '1x' 30 '1x' 30};
    % 定义列宽：固定宽度列 + 最后一列 '1x' 自动填充
    g.ColumnWidth = {100 100 100 100 100 100 100 100 '1x'};
    
    % 参数控件 (第一行)
    uilabel(g,'Text','Alpha (suppression)','HorizontalAlignment','right'); % 抑制强度 α
    sldAlpha = uislider(g,'Limits',[0 1],'Value',0.7); sldAlpha.Layout.Column=[2 4];
    
    uilabel(g,'Text','Mask smooth σ','HorizontalAlignment','right'); % 掩膜高斯平滑的标准差 sigma
    sldSigma = uislider(g,'Limits',[0 5],'Value',1.0); sldSigma.Layout.Column=[6 8];
    
    % 参数控件 (第二行)
    uilabel(g,'Text','Signal band (nm)','HorizontalAlignment','right'); % 用于 SNR 计算的信号波段
    edtSig = uieditfield(g,'text','Value','[868 905]'); edtSig.Layout.Column=[2 3];
    
    uilabel(g,'Text','Noise band (nm)','HorizontalAlignment','right'); % 用于 SNR 计算的噪声波段
    edtNoise = uieditfield(g,'text','Value','[905 inf]'); edtNoise.Layout.Column=[4 5];
    
    % 按钮控件 (第四行)
    btnBrush = uibutton(g,'Text','Freehand Brush','ButtonPushedFcn',@onBrush);  % 自由涂抹
    btnRect = uibutton(g,'Text','Add Rectangle','ButtonPushedFcn',@onRect); % 添加矩形
    btnPoly = uibutton(g,'Text','Add Polygon','ButtonPushedFcn',@onPoly); % 多边形（多顶点连成）
    btnUndo = uibutton(g,'Text','Undo','ButtonPushedFcn',@onUndo); % 撤销
    
    % 绘图工具按钮布局 (修正：将它们放在第 2 行，避免与滑块冲突)
    btnBrush.Layout.Row=2; btnBrush.Layout.Column=6; 
    btnRect.Layout.Row=2; btnRect.Layout.Column=7; 
    btnPoly.Layout.Row=2; btnPoly.Layout.Column=8; 
    btnUndo.Layout.Row=2; btnUndo.Layout.Column=9;
    
    btnClear = uibutton(g,'Text','Clear Mask','ButtonPushedFcn',@onClear);
    btnApply = uibutton(g,'Text','Apply Mask','ButtonPushedFcn',@onApply);
    btnExport = uibutton(g,'Text','Export Mask & Result','ButtonPushedFcn',@onExport);
    
    btnClear.Layout.Row=4; btnClear.Layout.Column=[1 3]; 
    btnApply.Layout.Row=4; btnApply.Layout.Column=[4 6]; 
    btnExport.Layout.Row=4; btnExport.Layout.Column=[7 9]; 
    
    % 轴区：绘图区域
    axFW = uiaxes(g); axFW.Layout.Row=3; axFW.Layout.Column=[1 5]; % F–W 热图 (CWT 幅度)
    axMask = uiaxes(g); axMask.Layout.Row=3; axMask.Layout.Column=[6 9]; % 掩膜预览
    axSpec = uiaxes(g); axSpec.Layout.Row=5; axSpec.Layout.Column=[1 9]; % 原始/重构光谱对比
    
    % 指标和提示 (第六行)
    lblSNR = uilabel(g,'Text','SNR: N/A'); lblSNR.Layout.Row=6; lblSNR.Layout.Column=[1 3]; % 信噪比显示
    lblHint= uilabel(g,'Text','Hint: Use Freehand/Rectangle/Polygon to add mask on F–W map, then "Apply".');
    lblHint.Layout.Row=6; lblHint.Layout.Column=[4 9]; % 操作提示
    
    % ---------------- 计算初值 ----------------
    % 连续小波变换 (CWT)
    [cfs,frq] = cwt(spectrum,'amor'); % cfs: CWT 系数 (F×N), frq: 归一化频率刻度 (F×1)
    % W 是 CWT 系数的归一化幅度，用于 F–W 热图显示
    W = abs(cfs); W = W ./ max(W(:)+eps);
    signalMean = mean(spectrum);
    
    % 掩膜初始化
    M = zeros(size(W),'logical'); % F×N 逻辑掩膜
    maskStack = {}; % 掩膜栈，用于实现撤销 (最多存 20 步)
    
    % 将关键数据存储到 App 对象的 userData 中，以确保回调函数可以访问
    ui.UserData.N = N;             % 数据点数量
    ui.UserData.cfs = cfs;         % 原始 CWT 系数
    ui.UserData.frq = frq;         % 频率刻度
    ui.UserData.W = W;             % 归一化 CWT 幅度
    ui.UserData.M = M;             % 当前二值掩膜
    ui.UserData.maskStack = maskStack; % 撤销栈
    ui.UserData.signalMean = signalMean;
    
    % 初始绘图：F–W 热图、掩膜预览、原始光谱对比（初始重构谱与原始谱相同）
    plotFW(); plotMask(); plotSpec(spectrum, spectrum);
    
    ui.UserData.sldAlpha = sldAlpha;
    ui.UserData.sldSigma = sldSigma;
    ui.UserData.edtSig = edtSig;
    ui.UserData.edtNoise = edtNoise;
    ui.UserData.lblSNR = lblSNR;
    
    % ---------------- 回调函数 ----------------
        function onBrush(~,~)
            % 在 FW 热图上自由涂抹，生成二值掩膜 (使用 drawfreehand)
            figure(ui); % 保持 App 窗口焦点
            roi = drawfreehand(axFW,'Color',[1 0 0],'LineWidth',1.5);
            addMaskFromROI(roi);
        end
    
        function onRect(~,~)
            % 在 FW 热图上绘制矩形，生成二值掩膜 (使用 drawrectangle)
            roi = drawrectangle(axFW,'Color',[0 0.6 1],'LineWidth',1.5);
            addMaskFromROI(roi);
        end
    
        function onPoly(~,~)
            % 绘制多边形，按点连接，生成二值掩膜
            figure(ui);
            % drawpolygon 允许用户通过点击来定义顶点。
            % 双击或右键点击最后一个点即可完成绘制并自动闭合。
            roi = drawpolygon(axFW, 'Color', [0.8 0.5 0], 'LineWidth', 1.5);
            addMaskFromROI(roi);
        end
    
        function addMaskFromROI(roi)
            M = ui.UserData.M;
            % 将 ROI (Region of Interest) 映射为 F×N 栅格上的逻辑掩膜
            BW = createMask(roi); % createMask 返回与 W 同大小的逻辑矩阵
            pushUndo(); % 压入当前 M 到撤销栈
            M = M | BW; % 逻辑或，将新区域添加到当前掩膜 M
            ui.UserData.M = M;
            plotMask(); % 更新掩膜预览图
        end
    
        function onUndo(~,~)
            % 撤销：弹出上一步的掩膜 M
            maskStack = ui.UserData.maskStack;
            
            if ~isempty(maskStack)
                % 2. 修改数据
                M = maskStack{end}; maskStack(end) = [];
                
                % 3. 存储回 UserData
                ui.UserData.M = M;
                ui.UserData.maskStack = maskStack;
                plotMask();
            end
        end
    
        function onClear(~,~)
            % 1. 获取共享数据
            M = ui.UserData.M;
            
            pushUndo(); % pushUndo 已经处理了 maskStack
            
            % 2. 修改数据
            M(:) = false;
            
            % 3. 存储回 UserData
            ui.UserData.M = M;
            plotMask();
        end
    
        function onApply(~,~)
            % 1. 获取共享数据
            cfs = ui.UserData.cfs;
           
            M = ui.UserData.M;
            signalMean = ui.UserData.signalMean;   % *** 周期范围 [min, max] ***
                
            % 2. 获取 UI 控件句柄和值
            sldAlpha_h = ui.UserData.sldAlpha;
            sldSigma_h = ui.UserData.sldSigma;
            lblSNR_h = ui.UserData.lblSNR;
        
            alpha = sldAlpha_h.Value; 
            sig = max(0.0, sldSigma_h.Value); 
        
            if sig>0
                Ms = imgaussfilt(double(M), sig);
                Ms = Ms./max(Ms(:)+eps);
            else
                Ms = double(M);
            end
            Ms = min(max(Ms,0),1);  % clamp
        
            cfs_new = cfs .* (1 - alpha*Ms);
            % === 5. ——频率归一化补偿—— ===
            % 对每个波长位置归一化，保持能量平衡
            w_orig = sum(abs(cfs), 1);
            w_new = sum(abs(cfs_new), 1);
            scale = w_orig ./ max(w_new, eps);
            scale = min(scale, 5);  % 防止放大过度
            cfs_eq = cfs_new .* scale;
            % *** 修正：向 icwt 传递频率刻度 frq ***
            % 格式: icwt(coefficients, wavelet_name, frequency_values, num_data_points)
            spec_rec = icwt(cfs_new, 'SignalMean', signalMean);
            spec_rec = spec_rec(:).';
            % === 7. 基线趋势校正 ===
            % 去掉掩膜引入的低频斜率，同时保持原信号均值
          
            r = spectrum - spec_rec;
            N = numel(r);
            % 估计要保留的DCT低阶数K
            K = max(3, floor(N/200));          % 极保守的下限(没有frq时)
            try
                % 如果 frq 可用且维度匹配，用掩膜中的最低频来自适应K
                if ~isempty(frq)
                    frq_vec = frq(:);
                    nfreq = size(Ms,1);
                    if numel(frq_vec) == nfreq
                        row_mask = any(Ms > 1e-6, 2);      % 哪些频率层被掩膜到
                        f_mask = frq_vec(row_mask);
                        if isempty(f_mask), f_mask = frq_vec; end
                        f_min = max(min(f_mask), 1e-9);    % 掩膜涉及的最低频
                        f_c   = 0.5 * f_min;               % 截止取其一半(更低)
                        % 经验换算：k ≈ 2*N*f (若frq已归一化到[0,0.5]则成立；
                        % 否则K会更小——更保守，只回填更慢的趋势)
                        K_adapt = max(1, floor(2 * N * f_c));
                        % 保护上限：只留极低频基线(不会回带条纹)
                        K = min(max(K, K_adapt), max(10, floor(N/150)));
                    end
                end
            catch
                % 忽略自适应失败，使用保守K
            end
        
            % 执行低阶DCT回填
            c = dct(r(:));
            c(K+1:end) = 0;
            trend_low = idct(c);
            trend_low = trend_low(:).';
        
            spec_rec = spec_rec + trend_low;                      % 只加回超低频基线
            spec_rec_corr = spec_rec + (signalMean - mean(spec_rec));  % 精确保持均值

            % ---- 线性幅度标定（保持均值不变）----
            if exist('lambda','var')
                refMask = (lambda >= 600) & (lambda <= 1800);   % 参考区间（你可换成 signalBand）
            else
                refMask = true(size(spec_rec_corr));
            end
            
            x_ref = spectrum(refMask);
            y_ref = spec_rec_corr(refMask);
            
            % 用分位数代替 min/max，更稳健
            x_lo = prctile(x_ref, 5);  x_hi = prctile(x_ref, 95);
            y_lo = prctile(y_ref, 5);  y_hi = prctile(y_ref, 95);
            
            % 斜率 a：对齐幅度；偏置 b：保证最终均值仍等于 signalMean
            a = (x_hi - x_lo) / max(y_hi - y_lo, eps);
            b = signalMean - a * mean(spec_rec_corr);
            
            spec_rec_corr = a * spec_rec_corr + b;  % 线性校正（均值保持不变）
            
            plotSpec(spectrum, spec_rec_corr);
            [snrVal, note] = bandSNR(lambda, spec_rec_corr);
        
            % 3. 更新 SNR 标签
            lblSNR_h.Text = sprintf('SNR (band-limited): %.2f dB  %s', snrVal, note);

            % === 9. 保存至 UserData ===
            ui.UserData.last_spec_rec = spec_rec_corr;
        end
    
        function onExport(~,~)
            % 1. 获取共享数据和控件句柄
            cfs = ui.UserData.cfs;
            M = ui.UserData.M;
            signalMean = ui.UserData.signalMean;
            frq = ui.UserData.frq; % *** 新增：获取频率刻度 ***
            sldAlpha_h = ui.UserData.sldAlpha;
            sldSigma_h = ui.UserData.sldSigma;
            edtSig_h = ui.UserData.edtSig;
            edtNoise_h = ui.UserData.edtNoise;
        
            [file,path] = uiputfile('mask_and_result.mat','Save Mask & Result');
            if isequal(file,0), return; end
        
            % 重新计算重构光谱（确保与当前参数一致）
            alpha = sldAlpha_h.Value; 
            sig = max(0.0, sldSigma_h.Value);
        
            % ===== 2) 掩膜平滑 =====
            if sig > 0
                Ms = imgaussfilt(double(M), sig);
                Ms = Ms ./ max(Ms(:) + eps);
            else
                Ms = double(M);
            end
            Ms = min(max(Ms,0),1);  % clamp
        
            % ===== 3) 核心三步(按你的写法) =====
            cfs_new  = cfs .* (1 - alpha*Ms);
            spec_rec = icwt(cfs_new, 'SignalMean', signalMean);
            spec_rec = spec_rec(:).';
        
            % ===== 4) 仅回填“超低频基线”(DCT低阶投影) =====
            r = spectrum - spec_rec;          % 残差(基线差)
            N = numel(r);
        
            % 基础(保守)低阶数
            K = max(3, floor(N/200));
        
            % 若 frq 与频率维匹配，则用掩膜中的最低频自适应K
            try
                if ~isempty(frq)
                    frq_vec = frq(:);
                    nfreq = size(Ms,1);
                    if numel(frq_vec) == nfreq
                        row_mask = any(Ms > 1e-6, 2);   % 哪些频率层被掩膜
                        f_mask = frq_vec(row_mask);
                        if isempty(f_mask), f_mask = frq_vec; end
                        f_min = max(min(f_mask), 1e-9);
                        f_c   = 0.5 * f_min;           % 截止频率=最低掩膜频率的一半
                        K_adapt = max(1, floor(2 * N * f_c));  % 经验换算: k≈2Nf
                        K = min(max(K, K_adapt), max(10, floor(N/150))); % 上限保护
                    end
                end
            catch
                % 忽略自适应失败，保守K照用
            end
        
            % 低阶DCT回填（只保留超低频）
            c = dct(r(:));
            c(K+1:end) = 0;
            trend_low = idct(c).';
            spec_rec = spec_rec + trend_low;
        
            % 均值严格对齐到 signalMean（数值消差）
            spec_rec = spec_rec + (signalMean - mean(spec_rec));

            % ===== 线性幅度标定（保持均值不变）=====
            % 参考区间：优先使用“Signal band”输入；否则退回 [600, 1800]；再否则全局
            refMask = true(size(spec_rec));
            try
                if ~isempty(edtSig_h) && ~isempty(edtSig_h.Value)
                    sigBand = eval(edtSig_h.Value);                 % 例如 [868 905]
                    if numel(sigBand)==2
                        refMask = (lambda >= sigBand(1)) & (lambda <= sigBand(2));
                    end
                elseif ~isempty(lambda)
                    refMask = (lambda >= 600) & (lambda <= 1800);
                end
            catch
                % 出错则使用默认（全局）
                refMask = true(size(spec_rec));
            end
        
            x_ref = spectrum(refMask);
            y_ref = spec_rec(refMask);
        
            % 用 5–95% 分位代替 min/max（避免尖峰干扰）
            x_lo = prctile(x_ref, 5);  x_hi = prctile(x_ref, 95);
            y_lo = prctile(y_ref, 5);  y_hi = prctile(y_ref, 95);
        
            a = (x_hi - x_lo) / max(y_hi - y_lo, eps);     % 斜率：幅度对齐
            b = signalMean - a * mean(spec_rec);           % 偏置：保持均值=signalMean
            spec_rec = a * spec_rec + b;                   % 线性校正（不改变均值）
        
            % 保存参数和数据
            params.alpha = alpha; params.sigma = sig;
            params.signalBand = eval(edtSig_h.Value);    % 通过句柄访问值
            params.noiseBand = eval(edtNoise_h.Value);  % 通过句柄访问值
        
            % 按是否存在 frq 一并保存
            if isempty(frq)
                save(fullfile(path,file), 'lambda','spectrum','cfs','M','params','spec_rec','Ms','signalMean');
            else
                save(fullfile(path,file), 'lambda','spectrum','cfs','frq','M','params','spec_rec','Ms','signalMean');
            end
        
            uialert(ui, 'Saved mask & reconstructed spectrum.', 'Export', 'Icon', 'success');
        end
    
    % ---------------- 工具函数 ----------------
        function plotFW()
            % 绘制 F–W 热图 (CWT 幅度)
            % 1. 获取共享数据
            W = ui.UserData.W;
            frq = ui.UserData.frq;
            
            imagesc(axFW, lambda, frq, W); axis(axFW,'xy'); 
            xlabel(axFW,'Wavelength (nm)'); ylabel(axFW,'Normalized frequency');
            title(axFW,'F–W spectra (|CWT|, normalized)'); colormap(axFW,jet); colorbar(axFW);

            axFW.XLim = [600 1800];          % 等价于 xlim(axSpec,[600 1800])
            axFW.XLimMode = 'manual';        % 防止后续自动重设（可选）

            axFW.YLim = [0 0.35];          % 等价于 xlim(axSpec,[600 1800])
            axFW.YLimMode = 'manual';        % 防止后续自动重设（可选）
        end
    
        function plotMask()
            % 绘制当前二值掩膜预览图
            % 1. 获取共享数据
            M = ui.UserData.M;
            frq = ui.UserData.frq;
            
            imagesc(axMask, lambda, frq, M); axis(axMask,'xy');
            xlabel(axMask,'Wavelength (nm)'); ylabel(axMask,'Normalized frequency');
            title(axMask,'Current mask (white=masked)');
            colormap(axMask,gray);
            axMask.XLim = [600 1800];          % 等价于 xlim(axSpec,[600 1800])
            axMask.XLimMode = 'manual';        % 防止后续自动重设（可选）
            axMask.YLim = [0 0.35];          % 等价于 xlim(axSpec,[600 1800])
            axMask.YLimMode = 'manual';        % 防止后续自动重设（可选）
        end
    
        function plotSpec(raw, rec)
            % 绘制原始光谱与重构光谱的对比图
            cla(axSpec); % 清空轴
            plot(axSpec, lambda, raw, 'Color',[0.5 0.5 0.5], 'LineWidth',1.7); hold(axSpec,'on');
            plot(axSpec, lambda, rec, 'b','LineWidth',2.0); hold(axSpec,'off');
            legend(axSpec, {'Raw','Masked-reconstructed'},'Location','best');
            xlabel(axSpec,'Wavelength (nm)'); ylabel(axSpec,'Intensity (a.u.)');
            title(axSpec,'Spectrum before/after masking');
            grid(axSpec,'on');
            axSpec.XLim = [600 1800];          % 等价于 xlim(axSpec,[600 1800])
            axSpec.XLimMode = 'manual';        % 防止后续自动重设（可选）
        end
    
        function pushUndo()
            % 压入当前掩膜到栈，并限制栈大小
            % 1. 获取共享数据
            M = ui.UserData.M;
            maskStack = ui.UserData.maskStack;
            
            % 2. 修改数据
            maskStack{end+1} = M;
            if numel(maskStack)>20, maskStack(1)=[]; end
            
            % 3. 存储回 UserData
            ui.UserData.maskStack = maskStack;
        end
    
        function [snrDB, note] = bandSNR(lam, spec)
            % 1. 获取控件句柄
            edtSig_h = ui.UserData.edtSig;
            edtNoise_h = ui.UserData.edtNoise;
        
            % 计算带限信噪比 (Band-limited Signal-to-Noise Ratio)
            try
                sigBand = eval(edtSig_h.Value); % 通过句柄访问值
                noiseBand = eval(edtNoise_h.Value); % 通过句柄访问值
            catch
                % 若用户输入无效，使用默认值
                sigBand = [868 905]; noiseBand = [905 inf];
            end
            
            % 光谱归一化到 [0, 1]
            s = (spec - min(spec)) / max(eps, (max(spec)-min(spec)));
            
            % 信号功率 (均方)
            idxSig = lam>=sigBand(1) & lam<=sigBand(2);
            P_signal = mean(s(idxSig).^2);
        
            % 噪声功率：在噪声带内，使用多项式基线下包络拟合残差作为噪声 (更精确的估计)
            idxNoise = lam>=noiseBand(1) & lam<=min(noiseBand(2), max(lam));
            x = find(idxNoise);
            y = s(idxNoise);
            
            if numel(x) > 20
                % 使用多次多项式拟合进行下包络拟合，以去除噪声带内可能存在的基线信号
                p = 6; xn = (x-mean(x))/std(x); yfit = y; % p=6 次多项式
                for k=1:100 % 迭代拟合，趋近于下包络
                    c = polyfit(xn, yfit, p);
                    yf = polyval(c, xn);
                    yfit = min(yfit, yf); % 关键步骤：取当前拟合和历史拟合的最小值 (下包络)
                end
                eta = y - yfit; % 噪声/残差 = 信号 - 基线/包络
                P_noise = mean(eta.^2); % 噪声功率 (均方残差)
                note = sprintf('[%.0f–%.0f] vs >%.0f nm',sigBand(1),sigBand(2),noiseBand(1));
            else
                % 数据点不足，直接使用噪声带内光谱的均方值作为噪声功率 (简化/粗略版)
                P_noise = mean(s(idxNoise).^2);
                note = sprintf('[%.0f–%.0f] vs >%.0f nm (raw)',sigBand(1),sigBand(2),noiseBand(1));
            end
            
            % SNR (dB) = 10 * log10 (P_signal / P_noise)
            snrDB = 10*log10(max(P_signal,eps)/max(P_noise,eps));
        end
end
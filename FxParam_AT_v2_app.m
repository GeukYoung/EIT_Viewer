function [Result,DataSet] = FxParam_AT_v2_app(DataSet)
% input
%   data : EIT data 
% output
%   tidal : tidal volume & time info

% DataSet.EIT.fs = round(1/seconds(median(diff(DataSet.EIT.t_hms))));
DataSet.EIT.fs = 100;
opt.th_lung = 0.2;
opt.n_avg = 3;
opt.sqi_dFRC = 0.3; % 0.2
opt.sqi_dRR = 0.5; % 0.2
opt.sqi_dTV = 2; % 1.5
% opt.sqi_TV = [10 2000];
opt.sqi_wCV = 100 * 0.2; % 200 ms
opt.sqi_wth = 100 * 60; % 10 s

% ESU
temp_ESU = movstd(DataSet.EIT.ohm',[opt.sqi_wCV,0])'./movmean(DataSet.EIT.ohm',[opt.sqi_wCV,0])';
temp_ESU = movmean(max(temp_ESU),[opt.sqi_wCV 0]);
th_sqi = movmin(temp_ESU,[opt.sqi_wth 0],'omitnan')*10;
sqi.ESU_100fs = temp_ESU<th_sqi;

if DataSet.EIT.Settings(end).Elec_config ~= 0
    switch DataSet.EIT.Settings(end).Elec_config
        case 0 % front-down
            num_elec = 1:16; % default
        case 1 % front-up
            num_elec = 16:-1:1;
        case 2 % back-down
            num_elec = [9:1:16 1:8];
        case 3 % back-up
            num_elec = [8:-1:1 16:-1:9];
        case 4 % side-left
            num_elec = [5:1:16 1:4];
        case 5 % side-right
            num_elec = [13:1:16 1:12];
    end
    [ro] = FxEIT_vReorder(num_elec);
    if ~isfield(DataSet.EIT,'ohm_HPF')
        DataSet.EIT.ohm = DataSet.EIT.ohm(ro,:);
    else
        DataSet.EIT.ohm = DataSet.EIT.ohm(ro,:);
        DataSet.EIT.ohm_HPF = DataSet.EIT.ohm_HPF(ro,:);
    end
end

%% HPF
if ~isfield(DataSet.EIT,'ohm_HPF')
    Wn = 0.7; N = 2;
    [b,a] = butter(N, Wn*2/DataSet.EIT.fs,'low');
    
    idx_notfin = ~isfinite(DataSet.EIT.ohm);
    if sum(sum(~idx_notfin,2)==0) > 0
        errordlg('Data is not availble (NaN)')
        return;
    else
%         DataSet.EIT.ohm_HPF = fillmissing(DataSet.EIT.ohm','nearest','EndValues','nearest')';
        DataSet.EIT.ohm_HPF = DataSet.EIT.ohm;
        DataSet.EIT.ohm_HPF(idx_notfin) = 0;
        DataSet.EIT.ohm_HPF = filtfilt(b,a,DataSet.EIT.ohm_HPF')';
        DataSet.EIT.ohm_HPF(idx_notfin) = NaN;
        DataSet.EIT.ohm = [];
    end
    clear a b Wn N idx_notfin;
end

%% RVS
if ~isfield(DataSet.EIT,'RVS')
    DataSet.EIT.RVS = FxRecon_AT(DataSet.EIT.ohm_HPF,'sum');
end

%% Check Saturation
opt.th_sat = 15;
opt.w_sat = 2*100;
DataSet.EIT.Sat = (mean(DataSet.EIT.ohm_HPF) > opt.th_sat) | isnan(mean(DataSet.EIT.ohm_HPF));
DataSet.EIT.Sat = movmax(DataSet.EIT.Sat,[opt.w_sat opt.w_sat]);
DataSet.EIT.RVS(DataSet.EIT.Sat) = nan;

%% Pressure index match
if sum(diff(DataSet.Pressure.idx_press) < 0)
    % find idx_EIT broken 
    tp_break_EIT = find(diff(DataSet.EIT.idx_scan) < 0); % find each end idx point
    n_break = DataSet.EIT.idx_scan(tp_break_EIT);
    n_break = cumsum(n_break);
    tp_break_EIT = [tp_break_EIT+1 length(DataSet.EIT.idx_scan)];
    
    % find idx_Pressure broken
    tp_break_P = find(diff(DataSet.Pressure.idx_press) < 0);
    tp_break_P = [tp_break_P+1 length(DataSet.Pressure.idx_press)];
    
    % find idx_ECG broken
    % AT(...?)
    for cnt = 1:length(tp_break_EIT)-1
        DataSet.Pressure.idx_press(tp_break_P(cnt):tp_break_P(cnt+1)) = DataSet.Pressure.idx_press(tp_break_P(cnt):tp_break_P(cnt+1))+n_break(cnt);
        DataSet.EIT.idx_scan(tp_break_EIT(cnt):tp_break_EIT(cnt+1)) = DataSet.EIT.idx_scan(tp_break_EIT(cnt):tp_break_EIT(cnt+1))+n_break(cnt);
    end
end

if DataSet.Pressure.idx_press(1) == 0
    DataSet.Pressure.idx_press = DataSet.Pressure.idx_press + 1;
end
DataSet.Pressure.Paw_origin = DataSet.Pressure.Paw;
DataSet.Pressure.Paw = zeros(size(DataSet.Pressure.idx_press));
DataSet.Pressure.Paw(DataSet.Pressure.idx_press(1:length(DataSet.Pressure.Paw_origin))) = DataSet.Pressure.Paw_origin;
DataSet.Pressure.Faw_origin = DataSet.Pressure.Faw;
DataSet.Pressure.Faw = zeros(size(DataSet.Pressure.idx_press));
DataSet.Pressure.Faw(DataSet.Pressure.idx_press(1:length(DataSet.Pressure.Faw_origin))) = DataSet.Pressure.Faw_origin;

%% Find peaks
peaks = sFxEIT_TVpeak(DataSet.EIT.RVS,DataSet.EIT.fs);
if isempty(peaks)
    Result = [];
    disp('peak does not found!');
    return;
end
DataSet.EIT.peaks = peaks;

for cnt = 1:length(peaks.insp)-1
    Result(cnt).idx_exp = peaks.exp(cnt);
    Result(cnt).idx_insp = peaks.insp(cnt);
    Result(cnt).t_exp = DataSet.EIT.t_hms(peaks.exp(cnt));
    Result(cnt).t_insp = DataSet.EIT.t_hms(peaks.insp(cnt));
    Result(cnt).td_exp = DataSet.EIT.t_hms(peaks.exp(cnt+1)) - DataSet.EIT.t_hms(peaks.insp(cnt));
    Result(cnt).td_insp = DataSet.EIT.t_hms(peaks.insp(cnt)) - DataSet.EIT.t_hms(peaks.exp(cnt));
    Result(cnt).td_tidal = Result(cnt).td_insp + Result(cnt).td_exp;
    Result(cnt).RR = 60/seconds(Result(cnt).td_tidal);
    Result(cnt).IE = seconds(Result(cnt).td_insp)/seconds(Result(cnt).td_tidal);
    Result(cnt).TVi_eit = DataSet.EIT.RVS(peaks.insp(cnt)) - DataSet.EIT.RVS(peaks.exp(cnt));
    Result(cnt).TVe_eit = DataSet.EIT.RVS(peaks.insp(cnt)) - DataSet.EIT.RVS(peaks.exp(cnt+1));
    Result(cnt).MV_eit = Result(cnt).TVi_eit * Result(cnt).RR / 1000;
    Result(cnt).dFRC = Result(cnt).TVi_eit - Result(cnt).TVe_eit;
    Result(cnt).AU2mL = 1;
    Result(cnt).unit_V = 'AU';
    Result(cnt).Epad_config = DataSet.EIT.Settings.Elec_config;
    
    if isfield(DataSet,'Pressure')
        Result(cnt).PIP = max(DataSet.Pressure.Paw(peaks.exp(cnt):peaks.exp(cnt+1)));
        Result(cnt).PEEP = min(DataSet.Pressure.Paw(peaks.exp(cnt):peaks.exp(cnt+1)));
        Result(cnt).dP = Result(cnt).PIP - Result(cnt).PEEP;
        
        flow_temp = DataSet.Pressure.Faw(peaks.exp(cnt):peaks.exp(cnt+1));
        if length(flow_temp) > 10
            Result(cnt).flow_offset = mean(flow_temp([1:5 end-4:end]));
        else
            Result(cnt).flow_offset = mean(flow_temp);
        end
        
        Result(cnt).TVi = sum(abs(flow_temp(flow_temp>0)))/60/100;
        Result(cnt).TVe = sum(abs(flow_temp(flow_temp<0)))/60/100;
        Result(cnt).MV = Result(cnt).TVi * Result(cnt).RR;
        Result(cnt).dFRC = Result(cnt).TVi - Result(cnt).TVe;
        Result(cnt).Cdyn = Result(cnt).TVi/Result(cnt).dP;
        Result(cnt).AU2mL = Result(cnt).TVi/Result(cnt).TVi_eit;        
        Result(cnt).unit_V = 'mL';
        clear flow_temp;
    end
    
    if isfield(DataSet,'path')
        Result(cnt).path = DataSet.path;
    end
    
    if isfield(DataSet,'idx_session')
        Result(cnt).idx_session = DataSet.idx_session;
    end
    
    if isfield(DataSet,'idx_file')
        Result(cnt).idx_file = DataSet.idx_file;
    end
end

for cnt = 1:length(peaks.insp)-1
    %% TV
    Result(cnt).R_exp = DataSet.EIT.ohm_HPF(:,peaks.exp(cnt));
    Result(cnt).R_insp = DataSet.EIT.ohm_HPF(:,peaks.insp(cnt));
    Result(cnt).Im_TV = FxRecon_AT(-(DataSet.EIT.ohm_HPF(:,peaks.insp(cnt))-DataSet.EIT.ohm_HPF(:,peaks.exp(cnt)))*median([Result.AU2mL]));
    Result(cnt).TVi_eit2ml_same = Result(cnt).TVi_eit*median([Result.AU2mL]);
    Result(cnt).Im_mask_lung = Result(cnt).Im_TV > max(Result(cnt).Im_TV) * opt.th_lung;
%     imagesc(Show_ImAT(double(Result(cnt).Im_mask_lung)))
%     imagesc(Show_ImAT(Result(cnt).Im_TV))

    %% Cdyn
    if isfield(DataSet,'Pressure')
        Result(cnt).Im_Cdyn = Result(cnt).Im_TV ./ Result(cnt).dP;
    end
    
    %% RVD
    tp = [Result(cnt).idx_exp Result(cnt).idx_insp];
    [Result(cnt).Im_RVD, Result(cnt).sdRVD] = sFxEIT_RVD(DataSet,tp);
    
    %% Pendelluft volume
    tp = [Result(cnt).idx_exp Result(cnt).idx_insp peaks.exp(cnt+1)]; % peaks.exp(cnt+1) for one more valley info
    [Result(cnt).Im_Pendelluft, Result(cnt).Vpendelluft, Result(cnt).Im_tPendelluft, Result(cnt).Vpdf_gain, Result(cnt).Vpdf_loss] = sFxEIT_Pendelluft(DataSet,tp);
    if 1
        Result(cnt).Im_Pendelluft = Result(cnt).Im_Pendelluft.*Result(cnt).AU2mL;
        Result(cnt).Vpendelluft = Result(cnt).Vpendelluft.*Result(cnt).AU2mL;
        Result(cnt).Vpdf_gain = Result(cnt).Vpdf_gain.*Result(cnt).AU2mL;
        Result(cnt).Vpdf_loss = Result(cnt).Vpdf_loss.*Result(cnt).AU2mL;
    end
    
    %% Popen
    if isfield(DataSet,'Pressure')
        tp = [Result(cnt).idx_exp Result(cnt).idx_insp];
        [Result(cnt).Im_Popen, Result(cnt).sdPopen] = sFxEIT_Popen(DataSet,tp);
    end
    
    %% GI
    if isfield(DataSet,'Im_TV')
        Result(cnt).GI =  sFxEIT_GI(Result(cnt).Im_TV, FxRecon_AT('isbnd'));
    else
        Result(cnt).GI =  sFxEIT_GI(FxRecon_AT(Result(cnt).R_exp-Result(cnt).R_insp), FxRecon_AT('isbnd'));
    end

    %% CoVx,y
    if isfield(DataSet,'Im_TV')
        [Result(cnt).CoVx, Result(cnt).CoVy] = sFxEIT_CoV(Result(cnt).Im_TV);
    else
        [Result(cnt).CoVx, Result(cnt).CoVy] = sFxEIT_CoV(FxRecon_AT(Result(cnt).R_exp-Result(cnt).R_insp));
    end

    %% avg filter
    if cnt > opt.n_avg-1
        Result(cnt).sdRVD = mean([Result(cnt-opt.n_avg+1:cnt).sdRVD]);
        Result(cnt).sdPopen = mean([Result(cnt-opt.n_avg+1:cnt).sdPopen]);
        Result(cnt).GI = mean([Result(cnt-opt.n_avg+1:cnt).GI]);
        Result(cnt).CoVx = mean([Result(cnt-opt.n_avg+1:cnt).CoVx]);
        Result(cnt).CoVy = mean([Result(cnt-opt.n_avg+1:cnt).CoVy]);
        Result(cnt).MV = mean([Result(cnt-opt.n_avg+1:cnt).MV]);
    end
end

%% SQI
% dFRC
sqi.dFRC = abs([Result.TVi_eit]-[Result.TVe_eit])./(([Result.TVi_eit]+[Result.TVe_eit])/2);
sqi.dFRC = ~movmax((sqi.dFRC > opt.sqi_dFRC),[1 1]); % good
% dRR
sqi.dRR = opt.sqi_dRR > movstd([Result.RR],[1 1])./movmean([Result.RR],[1 1]);
% TV
% sqi.TV = opt.sqi_TV(2) > [Result.TVi_eit2ml_same] & [Result.TVi_eit2ml_same] > opt.sqi_TV(1);
sqi.dTV = opt.sqi_dTV > movstd([Result.TVi],[1 1])./movmean([Result.TVi],[1 1]);

% ESU
for cnt = 1:length(Result)
    sqi.ESU(cnt) = sqi.ESU_100fs(Result(cnt).idx_insp);
end

sqi.good = movmin(sqi.dFRC & sqi.dTV & sqi.dRR & sqi.ESU, [1 1]);
DataSet.sqi = sqi;

% save
for cnt = 1:length(Result)
    Result(cnt).sqi_dFRC = sqi.dFRC(cnt);
    Result(cnt).sqi_dRR = sqi.dRR(cnt);
    Result(cnt).sqi_dTV = sqi.dTV(cnt);
%     Result(cnt).sqi_TV = sqi.TV(cnt);
    Result(cnt).sqi_ESU = sqi.ESU(cnt);
end

%% Sub functions
    function [CoVx, CoVy] = sFxEIT_CoV(Im_TV)
        Im_TV(Im_TV<0) = 0;
        total_sum = sum(Im_TV);
        Im_TV = reshape(Im_TV,sqrt(length(Im_TV)),sqrt(length(Im_TV)))'; % convert vector -> matrix
        
        % CoVx
        temp = 0;
        for i = 1:size(Im_TV,2)
            temp = temp + sum(Im_TV(:,i)) * (i/size(Im_TV,2));
        end
        CoVx = temp/total_sum;
        
        % CoVy
        temp = 0;
        for i = 1:size(Im_TV,1)
            temp = temp + sum(Im_TV(i,:)) * (i/size(Im_TV,1));
        end
        CoVy = temp/total_sum;
    end

    function [GI] = sFxEIT_GI(Im_TV, isbnd)
        Im_TV(Im_TV<0) = 0;
        Im_TV(isbnd) = [];
        GI = (sum(abs(Im_TV-median(Im_TV))))/sum(Im_TV); 
    end

    function [Im_RVD, sdRVD] = sFxEIT_RVD(DataSet,tp)
        opt.th_lung = 0.25;
        opt.th_rvd = 0.4;
        
        Im_RVD = zeros(FxRecon_AT('nelem'),1);
        temp_TV = FxRecon_AT(-(DataSet.EIT.ohm_HPF(:,tp(2))-DataSet.EIT.ohm_HPF(:,tp(1))));
        Im_mask_lung = temp_TV > (max(temp_TV)*opt.th_lung);
        idx_lung = find(Im_mask_lung);
        wave_sigma = FxRecon_AT(-(DataSet.EIT.ohm_HPF(:,tp(1):tp(2))-DataSet.EIT.ohm_HPF(:,tp(1))));
        wave_sigma = wave_sigma(idx_lung,:);
        wave_sigma_sum = sum(wave_sigma);
        t_global = find(wave_sigma_sum>(max(wave_sigma_sum)*opt.th_rvd),1,'first')/DataSet.EIT.fs*1000; % s -> ms;
        for cnt_RVD = 1:length(idx_lung)
            t_open(cnt_RVD) = find(wave_sigma(cnt_RVD,:)>(max(wave_sigma(cnt_RVD,:))*opt.th_rvd),1,'first')/DataSet.EIT.fs*1000; % s -> ms
            Im_RVD(idx_lung(cnt_RVD)) = t_open(cnt_RVD) - t_global;
        end
        
        if sum(Im_mask_lung) == 0
            sdRVD = NaN;
        else
            sdRVD = std(t_open);
        end
    end
        
    function [Im_Popen, sdPopen] = sFxEIT_Popen(DataSet,tp)
        opt.th_lung = 0.25;
        opt.th_popen = 0.1;
        
        Im_Popen = zeros(FxRecon_AT('nelem'),1);
        temp_TV = FxRecon_AT(-(DataSet.EIT.ohm_HPF(:,tp(2))-DataSet.EIT.ohm_HPF(:,tp(1))));
        Im_mask_lung = temp_TV > (max(temp_TV)*opt.th_lung);
        idx_lung = find(Im_mask_lung);
        wave_sigma = FxRecon_AT(-(DataSet.EIT.ohm_HPF(:,tp(1):tp(2))-DataSet.EIT.ohm_HPF(:,tp(1))));
        wave_sigma = wave_sigma(idx_lung,:);
        
        for cnt_RVD = 1:length(idx_lung)
            idx_open(cnt_RVD) = find(wave_sigma(cnt_RVD,:)>(max(wave_sigma(cnt_RVD,:))*opt.th_popen),1,'first')-1;
            Im_Popen(idx_lung(cnt_RVD)) = DataSet.Pressure.Paw(tp(1) + idx_open(cnt_RVD));
        end
        
        if sum(Im_mask_lung) == 0
            sdPopen = NaN;
        else
            sdPopen = std(Im_Popen(idx_lung));
        end
    end
          
    function [Im_Pendelluft, Vpendelluft, Im_tPendelluft, Vpdf_gain, Vpdf_loss] = sFxEIT_Pendelluft(DataSet,tp)
        opt.th_lung = 0.25;
        
        Im_tPendelluft = zeros(FxRecon_AT('nelem'),1);
        Im_Pendelluft = zeros(FxRecon_AT('nelem'),1);
        temp_TV = FxRecon_AT(-(DataSet.EIT.ohm_HPF(:,tp(2))-DataSet.EIT.ohm_HPF(:,tp(1))));
        Im_mask_lung = temp_TV > (nanmax(temp_TV)*opt.th_lung);
        idx_lung = find(Im_mask_lung); % lung ROI
%         wave_sigma = DataSet.EIT.Image.RM(idx_lung,:)*-(DataSet.EIT.ohm_HPF(:,tp(1):tp(3))-DataSet.EIT.ohm_HPF(:,tp(1)));
        wave_sigma = FxRecon_AT(-(DataSet.EIT.ohm_HPF(:,tp(1):tp(3))-DataSet.EIT.ohm_HPF(:,tp(1))));
        [~, idx_sExpG] = max(sum(wave_sigma)); % global peak time
        for cnt_Pen = 1:length(idx_lung)
            [~, idx_sExpR] = max(wave_sigma(cnt_Pen,:));
            Im_Pendelluft(idx_lung(cnt_Pen)) = wave_sigma(cnt_Pen,idx_sExpR) - wave_sigma(cnt_Pen,idx_sExpG);
            Im_tPendelluft(idx_lung(cnt_Pen)) = (idx_sExpR - idx_sExpG)/DataSet.EIT.fs;
            if Im_tPendelluft(idx_lung(cnt_Pen)) < 0 % Gas Lost (after global peak)
                Im_Pendelluft(idx_lung(cnt_Pen)) =  -Im_Pendelluft(idx_lung(cnt_Pen)); % sign code
            end
        end
        Vpendelluft = sum(abs(Im_Pendelluft(Im_mask_lung)));
        Vpdf_gain = sum(Im_Pendelluft(Im_Pendelluft>0 & Im_mask_lung));
        Vpdf_loss = sum(Im_Pendelluft(Im_Pendelluft<0 & Im_mask_lung));
        % positive: Gained, negative: Lost
    end

    function peaks = sFxEIT_TVpeak(RVS,fs)
        ws = round(1*fs);
        MAC = zeros(length(RVS),1);
        MAC(1:ws) = mean(RVS(1:ws));
        if MAC(ws) < RVS(ws)
            state_itc = 1;
        else
            state_itc = -1;
        end
        margin = 0.1;
        cnt_itc = ws;
        cnt_peak = 1;
        cnt_v = 1;
        i_itc = 0;
        roc = 0;
        
        for i = (ws+1):length(RVS)
            if i < 5*ws + 1
                MAC(i) = mean(RVS((i-ws):i)) - state_itc*std(RVS((i-ws):i))*margin;
            else
                MAC(i) = mean(RVS((i-ws):i)) - state_itc*std(RVS((i-5*ws):i))*margin;
            end
            
            if state_itc == 1
                if MAC(i) >= RVS(i)
                    [~,temp_roc] = max(RVS((i-cnt_itc):i));
                    roc(cnt_peak) = i - cnt_itc + temp_roc - 1;
                    i_itc(cnt_peak) = i;
                    cnt_peak = cnt_peak + 1;
                    cnt_itc = 0;
                    state_itc = -1;
                end
            elseif state_itc == -1
                if MAC(i) <= RVS(i)
                    [~,temp_roc] = min(RVS((i-cnt_itc):i));
                    roc(cnt_peak) = i - cnt_itc + temp_roc - 1;
                    i_itc(cnt_peak) = i;
                    cnt_peak = cnt_peak + 1;
                    cnt_itc = 0;
                    state_itc = 1;
                end
            end
            cnt_itc = cnt_itc + 1;
        end
        
        locs_Peak2 = roc;
        if length(locs_Peak2) > 5
            if RVS(locs_Peak2(1)) > RVS(locs_Peak2(2))
                locs_Peak2(1) = [];
            end
            if mod(length(locs_Peak2),2) == 1
                locs_Peak2(end) = [];
            end

            tidal_th = 0.3;
            idx_rm = [];
            for i = 2:(length(locs_Peak2)/2)
                pre_tidal = abs(RVS(locs_Peak2(2*(i-1)-1)) - RVS(locs_Peak2(2*(i-1))));
                cur_tidal = abs(RVS(locs_Peak2(2*i-1)) - RVS(locs_Peak2(2*i)));
                if (tidal_th*pre_tidal) > cur_tidal
                    idx_rm = [idx_rm i];
                end
            end

            if length(idx_rm) > 1
                locs_Peak2([(2*idx_rm-1) (2*idx_rm)]) = [];
            end
            peaks.exp = locs_Peak2(1:2:end)';
            peaks.insp = locs_Peak2(2:2:end)';
        else
            peaks = [];
        end
    end
end
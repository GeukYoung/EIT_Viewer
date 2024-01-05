function [Result,DataSet] = FxParam_ATb2b_app(DataSet)
% %% Test line
% [DataSet] = FxImport_EITb2b;

%% Cable direction conv
if isfield(DataSet,'ScaleFactor')
    if isfield(DataSet.ScaleFactor,'CableDirection')
        switch DataSet.ScaleFactor.CableDirection(end)
            case 0
                num_elec = 1:16; % default
            case 1
                num_elec = 16:-1:1;
            case 2
                num_elec = [9:1:16 1:8];
            case 3
                num_elec = [8:-1:1 16:-1:9];
            case 4
                num_elec = [5:1:16 1:4];
            case 5
                num_elec = [13:1:16 1:12];
        end
        [ro] = FxEIT_vReorder(num_elec);
        DataSet.EIT.Zmag_peak = DataSet.EIT.Zmag_peak(ro,:);
        DataSet.EIT.Zmag_valley = DataSet.EIT.Zmag_valley(ro,:);
        DataSet.EIT.Zphase_peak = DataSet.EIT.Zphase_peak(ro,:);
        DataSet.EIT.Zphase_valley = DataSet.EIT.Zphase_valley(ro,:);
%         DataSet.EIT.OV_peak = DataSet.EIT.OV_peak(ro,:);
%         DataSet.EIT.OV_valley = DataSet.EIT.OV_valley(ro,:);
    end
end

%% mask 256->208
% DataSet.EIT.Zmag_peak(FxEIT_mask(16),:) = [];
% DataSet.EIT.Zmag_valley(FxEIT_mask(16),:) = [];
% DataSet.EIT.Zphase_peak(FxEIT_mask(16),:) = [];
% DataSet.EIT.Zphase_valley(FxEIT_mask(16),:) = [];
% DataSet.EIT.OV_peak(FxEIT_mask(16),:) = [];
% DataSet.EIT.OV_valley(FxEIT_mask(16),:) = [];

%% EIT FPV match
%   peak(n)       peak(n+1)  
%  ~~~^~~~~~~v~~~~~~^~~~
%         valley(n)
%  Breath(n) = valley(n)->peak(n+1) // need 1 peak delay for valley->peak (current peak->valley)
%  TVi(n) = TVi_flow(n+1)

idx_PFV = findidx(DataSet.PFV.TimeStamp,DataSet.EIT.TimeStamp);
idx_RVSp = findidx(DataSet.Wave.TimeStamp,DataSet.EIT.tRVSpeak);
idx_RVSv = findidx(DataSet.Wave.TimeStamp,DataSet.EIT.tRVSvalley);

% figure; plot(DataSet.Wave.TimeStamp,DataSet.Wave.RVS);
% hold on; plot(DataSet.EIT.tRVSpeak,ones(size(DataSet.EIT.tRVSpeak)),'o'); % missmatch

% scale factor setup
ScaleFactor_TV = repmat(double(DataSet.ScaleFactor.ScaleFactor(1)),size(DataSet.EIT.TimeStamp'));
if length(DataSet.ScaleFactor.ScaleFactor) > 1
    for cnt = 2:length(DataSet.ScaleFactor.ScaleFactor)
        ScaleFactor_TV(DataSet.EIT.TimeStamp>DataSet.ScaleFactor.TimeStamp(cnt)) = DataSet.ScaleFactor.ScaleFactor(cnt);
    end
end

%%
for cnt = 1:length(DataSet.EIT.TimeStamp)-1
    %% time info
    Result(cnt).t_exp = DataSet.EIT.tRVSvalley(cnt);
    Result(cnt).t_insp = DataSet.EIT.tRVSpeak(cnt+1);
    Result(cnt).td_insp = DataSet.EIT.tRVSpeak(cnt+1) - DataSet.EIT.tRVSvalley(cnt);
    Result(cnt).td_exp = DataSet.EIT.tRVSvalley(cnt+1) - DataSet.EIT.tRVSpeak(cnt+1);
    Result(cnt).td_tidal = Result(cnt).td_insp + Result(cnt).td_exp;
    Result(cnt).RR = 60/seconds(Result(cnt).td_tidal);
    Result(cnt).IE = Result(cnt).td_insp/Result(cnt).td_tidal;
    
    %% Volume EIT
    Result(cnt).ScaleFactor_TV = ScaleFactor_TV(cnt);
    Result(cnt).TVi_eit = Result(cnt).ScaleFactor_TV*FxRecon_AT((DataSet.EIT.Zmag_peak(:,cnt+1)-DataSet.EIT.Zmag_valley(:,cnt)),'sum')*1000; % mL
    Result(cnt).TVe_eit = Result(cnt).ScaleFactor_TV*FxRecon_AT((DataSet.EIT.Zmag_peak(:,cnt+1)-DataSet.EIT.Zmag_valley(:,cnt+1)),'sum')*1000; % mL
    Result(cnt).MV_eit = Result(cnt).TVi_eit * Result(cnt).RR / 1000; % L/cmH2O
    Result(cnt).dFRC_eit = Result(cnt).TVi_eit - Result(cnt).TVe_eit;

    %% recon image
    Result(cnt).Im_TV = FxRecon_AT(-(DataSet.EIT.Zmag_peak(:,cnt+1)-DataSet.EIT.Zmag_valley(:,cnt)));
%     Result(cnt).Im_TV = Result(cnt).ScaleFactor_TV*FxRecon_AT(-(DataSet.EIT.Zmag_peak(:,cnt+1)-DataSet.EIT.Zmag_valley(:,cnt)));
%     Show_ImAT(Result(cnt).Im_TV);

    %% flow sensor
    Result(cnt).PIP = DataSet.PFV.Ppeak(idx_PFV(cnt)); % cmH2O
    Result(cnt).PEEP = DataSet.PFV.PEEP(idx_PFV(cnt)); % cmH2O
    Result(cnt).dP = Result(cnt).PIP - Result(cnt).PEEP; % cmH2O
    Result(cnt).TVi = DataSet.PFV.TVi(idx_PFV(cnt)); % mL
    Result(cnt).TVe = DataSet.PFV.TVe(idx_PFV(cnt)); % mL
    Result(cnt).MV = Result(cnt).TVi * Result(cnt).RR/1000; % L
    Result(cnt).dFRC = Result(cnt).TVi - Result(cnt).TVe; % mL
    Result(cnt).Cdyn = Result(cnt).TVi/Result(cnt).dP; % mL/cmH2O
%     Result(cnt).AU2mL = Result(cnt).TVi/Result(cnt).TVi_au; 
%     Result(cnt).Epad_config = DataSet.EIT.Settings.Elec_config;
    
    %% FEIT
    % RVD
    Result(cnt).Im_RVD = DataSet.EIT.Im_RVD(:,cnt);
    Result(cnt).sdRVD = sFxEIT_sdRVD(Result(cnt).Im_RVD,Result(cnt).Im_TV);
    
    % GI
    Result(cnt).GI =  sFxEIT_GI(Result(cnt).Im_TV, FxRecon_AT('isbnd'));

    % CoVx,y
    [Result(cnt).CoVx, Result(cnt).CoVy] = sFxEIT_CoV(Result(cnt).Im_TV);
    
    % Im_TV AU -> V scale
    Result(cnt).Im_TV = Result(cnt).ScaleFactor_TV*Result(cnt).Im_TV;
    
    %% SQI
    Result(cnt).SQI = DataSet.Wave.SQI(idx_RVSp(cnt));
end

%% SQI
opt.th_lung = 0.2;
opt.n_avg = 3;
opt.sqi_dFRC = 0.3; % 0.2
opt.sqi_dRR = 0.5; % 0.2
opt.sqi_dTV = 2; % 1.5
% opt.sqi_TV = [10 2000];
opt.sqi_wCV = 100 * 0.2; % 200 ms
opt.sqi_wth = 100 * 10; % 10 s

% dFRC
sqi.dFRC = abs([Result.TVi_eit]-[Result.TVe_eit])./(([Result.TVi_eit]+[Result.TVe_eit])/2);
sqi.dFRC = ~movmax((sqi.dFRC > opt.sqi_dFRC),[1 1]); % good

% dRR
sqi.dRR = opt.sqi_dRR > movstd([Result.RR],[1 1])./movmean([Result.RR],[1 1]);

% TV
% sqi.TV = opt.sqi_TV(2) > [Result.TVi_au2ml_same] & [Result.TVi_au2ml_same] > opt.sqi_TV(1);
sqi.dTV = opt.sqi_dTV > movstd([Result.TVi],[1 1])./movmean([Result.TVi],[1 1]);

% % ESU
% temp_ESU = movstd(DataSet.EIT.ohm',[opt.sqi_wCV,0])'./movmean(DataSet.EIT.ohm',[opt.sqi_wCV,0])';
% temp_ESU = movmean(max(temp_ESU),[opt.sqi_wCV 0]);
% th_sqi = movmin(temp_ESU,[opt.sqi_wth 0],'omitnan')*10;
% sqi.ESU_100fs = temp_ESU<th_sqi;
% for cnt = 1:length(Result)
%     sqi.ESU(cnt) = sqi.ESU_100fs(Result(cnt).idx_insp);
% end

% save
for cnt = 1:length(Result)
    Result(cnt).sqi_dFRC = sqi.dFRC(cnt);
    Result(cnt).sqi_dRR = sqi.dRR(cnt);
    Result(cnt).sqi_dTV = sqi.dTV(cnt);
%     Result(cnt).sqi_TV = sqi.TV(cnt);
%     Result(cnt).sqi_ESU = sqi.ESU(cnt);
end

%% avg filter
% opt.n_avg = 5;
% if cnt > opt.n_avg-1
%     Result(cnt).sdRVD = mean([Result(cnt-opt.n_avg+1:cnt).sdRVD]);
%     Result(cnt).sdPopen = mean([Result(cnt-opt.n_avg+1:cnt).sdPopen]);
%     Result(cnt).GI = mean([Result(cnt-opt.n_avg+1:cnt).GI]);
%     Result(cnt).CoVx = mean([Result(cnt-opt.n_avg+1:cnt).CoVx]);
%     Result(cnt).CoVy = mean([Result(cnt-opt.n_avg+1:cnt).CoVy]);
%     Result(cnt).MV = mean([Result(cnt-opt.n_avg+1:cnt).MV]);
% end

%% Sub functions
    function [CoVx, CoVy] = sFxEIT_CoV(Im_TV)
        Im_TV(Im_TV<0) = 0;
        total_sum = sum(Im_TV);
        Im_TV = reshape(Im_TV,sqrt(length(Im_TV)),sqrt(length(Im_TV)))';
        
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

    function [sdRVD] = sFxEIT_sdRVD(Im_RVD,Im_TV)
        opt.th_lung = 0.25;
%         opt.th_rvd = 0.4;
        
        Im_mask_lung = Im_TV > (max(Im_TV)*opt.th_lung);
        idx_lung = Im_mask_lung==1;

        if sum(Im_mask_lung) == 0
            sdRVD = NaN;
        else
            sdRVD = std(Im_RVD(idx_lung),'omitnan');
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
            
    function [idx,lack] = findidx(A,tp)
        if isdatetime(A(1))
            A = datenum(A);
        end
        if isdatetime(tp(1))
            tp = datenum(tp);
        end
        idx = NaN(size(tp));
        lack = NaN(size(tp));
        cnt_2 = 1;
        for cnt_1 = 2:length(A)
            if A(cnt_1) > tp(cnt_2)
                if A(cnt_1)-tp(cnt_2) > A(cnt_1-1)-tp(cnt_2)
                    idx(cnt_2) = cnt_1-1;
                    lack(cnt_2) = A(cnt_1-1)-tp(cnt_2);
                else
                    idx(cnt_2) = cnt_1;
                    lack(cnt_2) = A(cnt_1)-tp(cnt_2);
                end
                cnt_2 = cnt_2 + 1;
            end
            if cnt_2 > length(tp)
                idx(isnan(idx)) = length(A);
                break;
            end
        end
        idx(isnan(idx)) = length(A);
    end
%     figure; plot(A,'o'); hold on;
%     plot(idx, A(idx), 'ro','MarkerFaceColor','r');
end
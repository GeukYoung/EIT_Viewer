function [Result,DataSet] = FxParam_HVb2b_app(DataSet)

figure; plot(DataSet.SV.TimeStamp,DataSet.SV.SV);

%% CVS
figure;
plot(DataSet.Wave.TimeStamp, DataSet.Wave.CVS); hold on;
plot(DataSet.peak_CVS.TimeStamp,zeros(size(DataSet.peak_CVS.TimeStamp)),'.');

%%
figure;
% h(1) = subplot(211);
plot(DataSet.Wave.TimeStamp, DataSet.Wave.CVS); hold on;
plot(DataSet.SV.TimeStamp,zeros(size(DataSet.SV.SV)),'.');
plot(DataSet.peak_CVS.tPeak_CVS,zeros(size(DataSet.peak_CVS.tPeak_CVS)),'rv');

% h(2) = subplot(212);
plot(DataSet.Wave.TimeStamp, movmax(DataSet.Wave.SQI,[10 0])); hold on;
linkaxes(h,'x');

figure;
plot(DataSet.Wave.TimeStamp, DataSet.Wave.ECG(1:5:end)); hold on;
plot(DataSet.SV.TimeStamp,zeros(size(DataSet.SV.SV)),'.');
plot(DataSet.peak_ECG.tPeak_ECG,zeros(size(DataSet.peak_ECG.tPeak_ECG)),'rv');


%%
idx_NIBP = findidx(DataSet.SV.TimeStamp,DataSet.NIBP.TimeStamp);
Result.SBP(idx_NIBP) = DataSet.NIBP.SBP;
Result.DBP(idx_NIBP) = DataSet.NIBP.DBP;
Result.MAP(idx_NIBP) = DataSet.NIBP.MAP;
Result.PR_NIBP(idx_NIBP) = DataSet.NIBP.PR_NIBP;

idx_Z = findidx(DataSet.SV.TimeStamp,DataSet.Z.TimeStamp);
Result.Z(idx_Z) = DataSet.Z.Z;

idx_SpO2 = findidx(DataSet.SV.TimeStamp,DataSet.SpO2.TimeStamp);
Result.SpO2(idx_SpO2) = DataSet.SpO2.SpO2;
Result.PR_SpO2(idx_SpO2) = DataSet.SpO2.PR_SpO2;



% scale factor setup
ScaleFactor_SV = repmat(double(DataSet.ScaleFactor.ScaleFactor(1)),size(DataSet.EIT.TimeStamp'));
if length(DataSet.ScaleFactor.ScaleFactor) > 1
    for cnt = 2:length(DataSet.ScaleFactor.ScaleFactor)
        ScaleFactor_SV(DataSet.EIT.TimeStamp>DataSet.ScaleFactor.TimeStamp(cnt)) = DataSet.ScaleFactor.ScaleFactor(cnt);
    end
end

%%
for cnt = 1:length(DataSet.SV.SV)-1
    %% time info
    Result(cnt).t_exp = DataSet.EIT.tRVSvalley(cnt);
    Result(cnt).t_insp = DataSet.EIT.tRVSpeak(cnt+1);
    Result(cnt).td_insp = DataSet.EIT.tRVSpeak(cnt+1) - DataSet.EIT.tRVSvalley(cnt);
    Result(cnt).td_exp = DataSet.EIT.tRVSvalley(cnt+1) - DataSet.EIT.tRVSpeak(cnt+1);
    Result(cnt).td_tidal = Result(cnt).td_insp + Result(cnt).td_exp;
    Result(cnt).RR = 60/seconds(Result(cnt).td_tidal);
    Result(cnt).IE = Result(cnt).td_insp/Result(cnt).td_tidal;
    
    %% Volume EIT
    Result(cnt).ScaleFactor_TV = ScaleFactor_SV(cnt);
    Result(cnt).TVi_eit = Result(cnt).ScaleFactor_TV*FxRecon_AT((DataSet.EIT.Zmag_peak(:,cnt+1)-DataSet.EIT.Zmag_valley(:,cnt)),'sum')*1000; % mL
    Result(cnt).TVe_eit = Result(cnt).ScaleFactor_TV*FxRecon_AT((DataSet.EIT.Zmag_peak(:,cnt+1)-DataSet.EIT.Zmag_valley(:,cnt+1)),'sum')*1000; % mL
    Result(cnt).MV_eit = Result(cnt).TVi_eit * Result(cnt).RR / 1000; % L/cmH2O
    Result(cnt).dFRC_eit = Result(cnt).TVi_eit - Result(cnt).TVe_eit;

    %% recon image
    Result(cnt).Im_TV = Result(cnt).ScaleFactor_TV*FxRecon_AT(-(DataSet.EIT.Zmag_peak(:,cnt+1)-DataSet.EIT.Zmag_valley(:,cnt)));
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
    function SVV = calcSVV
        temp;
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
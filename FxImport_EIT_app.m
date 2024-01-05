function [Data] = FxImport_EIT_app(varargin)
% - inputs
% filepath: '...\RawData'
% block: 'first'/'last'/[2 3 4]
% mask: 'inable'/'disable'/[1 2 16 17...]
% type: 'mag'/'comp'

if nargin < 1
    varargin{1} = uigetdir; 
end

if ischar(varargin{1})
    Path_main = varargin{1};
    if ~contains(Path_main, [filesep 'RawData'])
        temp_dir = dir(Path_main);
        if contains([temp_dir.name], 'RawData')
            Path_main = [Path_main [filesep 'RawData']];
        end
        clear temp_dir;
    end
    varargin = varargin(2:end);
%     argOffset = 1;
else
    error('ERROR: Data path');
end

opt.block = 0;
opt.mask = 1;
opt.type = 'magnitude';
opt.disp = 1;

numOrigInputArgs = numel(varargin);
if numOrigInputArgs ~= 0
    while numOrigInputArgs
        switch varargin{1}
            case 'block'
                if isnumeric(varargin{2})
                    opt.block = varargin{2};
                else
                    if contains(varargin{2},'first')
                        opt.block = 1;
                    elseif contains(varargin{2},'last')
                        opt.block = 99;
                    end
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'disp'
                if isnumeric(varargin{2})
                    opt.disp = varargin{2};
                else
                    error('disp');
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'mask'
                if contains(varargin{2},'disable')
                    opt.mask = 0;
                elseif contains(varargin{2},'inable')
                    opt.mask = 1;
                else
                    opt.mask = 1;
                    Data.ch_mask = varargin{2};
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'type'
                if contains(varargin{2},'mag') || contains(varargin{2},'abs')
                    opt.type = 'magnitude';
                elseif contains(varargin{2},'com')
                    opt.type = 'complex';
                elseif contains(varargin{2},'phase') || contains(varargin{2},'angle')
                    opt.type = 'phase';
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            otherwise
                msgbox(['Unvalid input :' varargin{1}]);
                numOrigInputArgs = 0;
        end
    end
end

if opt.disp
    f = waitbar(0,'Please Wait...','Name','Loading...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(f,'canceling',0);
end

dirlist = dir(Path_main);

cnt.eit = 0;
cnt.voltage = 0;
cnt.phase = 0;
cnt.ecg = 0;
cnt.cm = 0;
cnt.setting = 0;
cnt.pressure = 0;
cnt.pleth = 0;
cnt.spO2 = 0;

for cnt_file = 1:length(dirlist)
    temp_string = strsplit(dirlist(cnt_file).name,{'.', '['});
    if contains(temp_string{1},'EIT_','IgnoreCase',true)
        idx.eit(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.eit = cnt.eit + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).eit = fullfile(Path_main, dirlist(cnt_file).name);
        opt.flag_Z = 0;
    end
    
    if contains(temp_string{1},'EITImpedance_','IgnoreCase',true)
        idx.eit(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.eit = cnt.eit + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).eit = fullfile(Path_main, dirlist(cnt_file).name);
        opt.flag_Z = 1;
    end
    
    if contains(temp_string{1},'EITGain_','IgnoreCase',true)
        idx.eit(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.eit = cnt.eit + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).eit = fullfile(Path_main, dirlist(cnt_file).name);
        opt.flag_Z = 2;
    end
    
    if contains(temp_string{1},'Voltage','IgnoreCase',true)
        idx.voltage(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.voltage = cnt.voltage + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).voltage = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Phase','IgnoreCase',true)
        idx.phase(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.phase = cnt.phase + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).phase = fullfile(Path_main,dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'ECG','IgnoreCase',true)
        idx.ecg(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.ecg = cnt.ecg + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).ecg = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'EITCM','IgnoreCase',true)
        idx.cm(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.cm = cnt.cm + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).cm = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Setting','IgnoreCase',true)
        idx.setting(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.setting = cnt.setting + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).setting = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Pressure','IgnoreCase',true)
        idx.pressure(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.pressure = cnt.pressure + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).pressure = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Pleth','IgnoreCase',true)
        idx.pleth(str2double(regexp(temp_string{1},'\d*','match'))+1) = true;
        cnt.pleth = cnt.pleth + 1;
        Path(str2double(regexp(temp_string{1},'\d*','match'))+1).pleth = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'SpO2','IgnoreCase',true)
        idx_spO2 = regexp(temp_string{1},'\d*','match'); % cuz SpO[2] have number among name
        idx_spO2 = str2double(idx_spO2(end));
        idx.spO2(idx_spO2+1) = true;
        cnt.spO2 = cnt.spO2 + 1;
        Path(idx_spO2+1).spO2 = fullfile(Path_main, dirlist(cnt_file).name);
        clear idx_spO2;
    end 
end

total_EIT = cnt.eit; clear cnt_file;

if cnt.eit > 0
    if sum(idx.eit) ~= length(idx.eit)
        msgbox(['EIT_' num2str(find(idx.eit==0)) 'bin data missing']);
    end
end
if cnt.voltage > 0
    if sum(idx.voltage) ~= length(idx.voltage)
        msgbox(['Voltage_' num2str(find(idx.voltage==0)) 'data missing']);
    end
end
if cnt.phase > 0
    if sum(idx.phase) ~= length(idx.phase)  
        msgbox(['Phase_' num2str(find(idx.phase==0)) 'data missing']);
    end
end
if cnt.ecg > 0
    if sum(idx.ecg) ~= length(idx.ecg)
        msgbox(['ECG_' num2str(find(idx.ecg==0)) 'data missing']);
    end
end
if cnt.cm > 0
    if sum(idx.cm) ~= length(idx.cm)
        msgbox(['EITCM_ ' num2str(find(idx.cm==0)) 'data missing']);
    end
end
if cnt.pressure > 0
    if sum(idx.pressure) ~= length(idx.pressure)
        msgbox(['Pressure_ ' num2str(find(idx.pressure==0)) 'data missing']);
    end
end
if cnt.pleth > 0
    if sum(idx.pleth) ~= length(idx.pleth)
        msgbox(['Pleth_ ' num2str(find(idx.pleth==0)) 'data missing']);
    end
end
if cnt.spO2 > 0
    if sum(idx.spO2) ~= length(idx.spO2)
        msgbox(['SpO2_ ' num2str(find(idx.spO2==0)) 'data missing']);
    end
end

if opt.block ~= 0
    if size(opt.block) == 1
        if isfield(Path,'setting')
            Path_setting = Path(1).setting;
        end
        switch opt.block
            case 1
                Path = Path(1);
                total_EIT = 1;
            case 99
                Path = Path(end);
                total_EIT = 1;
            otherwise
                Path = Path(opt.block);
                total_EIT = 1;
        end
        if isfield(Path,'setting')
            Path(1).setting = Path_setting;
        end
    elseif length(opt.block) > 1
        try
            if isfield(Path,'setting')
                Path_setting = Path(1).setting;
            end
            Path = Path(opt.block);
            if isfield(Path,'setting')
                Path(1).setting = Path_setting;
            end
            total_EIT = length(opt.block);
        catch
            msgbox('block data is not match');
        end
    end
end

try
    % check data version
    fid = fopen(Path(1).eit);
    version = fread(fid,4,'uint8');
    version = version(4)*2^3 + version(3)*2^2 + version(2)*2^1 + version(1)*2^0;
    fclose(fid);
catch
    msgbox('Import failed');
%     fclose(fid);
end


Data.version = version;
Data.path = Path_main;

% AT HV info
if cnt.pressure > 0
    Data.SW_type = 'AT';
    Data.ch_mask = FxEIT_mask(16);
else
    Data.SW_type = 'HV';
    Data.ch_mask = [1;2;16;17;18;19;34;35;36;37;52;53;54;69;70;71;86;87;88;103;104;105;120;121;122;137;138;139;154;155;156;171;172;173;188;189;190;205;206;207;208;209;223;224;225;226;232;233;241;249;250;256];
end

% patient info
try
    fid = fopen(strcat(erase(Path_main,'\RawData'),'\Info\Patientinfo.ini'));
    info_p = textscan(fid,'%s');
    info_p = info_p{1,1};
    fclose(fid);

    for i = 1:length(info_p)
        info_name = split(info_p{i} ,'=');
        switch info_name{1}
            case 'ID'
                Data.info_p.ID = info_p{i}(length('ID')+2:end);
            case 'AGE'
                Data.info_p.Age = str2num(info_p{i}(length('AGE')+2:end));
            case 'GENDER'
                Data.info_p.Gender = info_p{i}(length('GENDER')+2:end);
            case 'WEIGHT'
                Data.info_p.W = str2num(info_p{i}(length('WEIGHT')+2:end));
            case 'HEIGHT'
                Data.info_p.H = str2num(info_p{i}(length('HEIGHT')+2:end));
        end
    end
catch
end

try
switch version
%% v0101 for AirTom/HemoVista
    case 5 % 0101      
        % raw data
        header.meas_num = 256; % 256 measurement
        header.header_num = 24; % 24 byte header of 1 data
        header.fs_EIT = 100;
        header.scan_num = 360000;
        header.protocolkey1 = 5;
        header.protocolkey2 = 6;
        header.ci_flag = 7;
        header.as_flag = 8;
        header.idx_scan = 9:12;
        header.timestamp = 13+1; % bug fixed(x) 13 + 1
        header.bValid = 23;
        header.bSQI = 24;
        
        fid = fopen(Path(1).eit);
        data_version = fread(fid,10000,'uint8');
        data_version = typecast(uint8(data_version),'int32');
        fclose(fid);
        
        % check data format
        if sum(diff(data_version(1:2328/4:end))) == 0 % Impedance (4 byte for 1 data)
            opt.data_ver = 1;
            header.data_num = 2328 - header.header_num;
        elseif sum(diff(data_version(1:2584/4:end))) == 0 % w/ gain
            opt.data_ver = 2;
            header.data_num = 2584 - header.header_num;
        elseif sum(diff(data_version(1:1536/4:end))) == 0 % pre
            opt.data_ver = 0;
            header.data_num = 1536 - header.header_num;
        else
            opt.data_ver = 0;
            header.data_num = 1536;
        end
        
        % setting file
        header.set_scale_volt = 0.0000221;
        header.set_total_num = 2632;
        header.set_size_interch = 2608;
        header.set_Amp = 1:8;
        header.set_Gain_table_CodeKey = 9;
        header.set_ScanIdx = 2618:2621;
        header.set_TimeStamp = (2622:2629)+1; % idx error +1
        header.set_GainX = 1:128;
        header.set_GainDigi = 129:144;
        header.set_GainDigi2 = 145:160;
        header.set_Protocoltype = 161;
        header.set_ProjNum = 162;
        header.set_Update = 163;
        
        % ecg file
        header.ecg_size = 69;
        header.ecg_idx = 2;
        header.ecg_data = 6;

        % pressure file
        header.pres_size = 16;
%         header.pre_idx = 2;
%         header.pre_data = 6;

        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            if opt.disp
                if getappdata(f,'canceling')
                    break;
                end
            end
            % EIT raw
            if 1
                if opt.disp
                    waitbar(cnt_EIT/total_EIT, f, ['Loading EIT bin data file '  num2str(cnt_EIT) ' of ' num2str(total_EIT)]);
                end
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            
            if rem(size(raw_data,1)/(header.header_num + header.data_num),10000) == 0
                header.scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            end
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            if mod(size(raw_data,1),(header.header_num + header.data_num)) ~= 0
                raw_data(floor(scan_num)*(header.header_num + header.data_num)+1:end) = [];
                scan_num = floor(scan_num);
            end
            if cnt_EIT < total_EIT && scan_num ~= header.scan_num
                msgbox(['data missing: EIT_bin' num2str(cnt_EIT) ' (' num2str(header.scan_num - scan_num) ' scan)']);
            end
            raw_data = reshape(raw_data,(header.header_num + header.data_num),scan_num);

%             header_data=raw_data(1:header.header_num,:);
            
            Data.EIT.ProtocolKey1(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.protocolkey1,:);
            Data.EIT.ProtocolKey2(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.protocolkey2,:);
            Data.EIT.CI_flag(cnt_scan:cnt_scan+scan_num-1) = logical(raw_data(header.ci_flag,:));
            Data.EIT.AS_flag(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.as_flag,:);   
            Data.EIT.idx_scan(cnt_scan:cnt_scan+scan_num-1) = raw_data(header.idx_scan(4),:)*256^3 + raw_data(header.idx_scan(3),:)*256^2 + raw_data(header.idx_scan(2),:)*256^1 + raw_data(header.idx_scan(1),:);
            Data.EIT.Valid(cnt_scan:cnt_scan+scan_num-1) = logical(raw_data(header.bValid,:));
            Data.EIT.SQI(cnt_scan:cnt_scan+scan_num-1) = logical(raw_data(header.bSQI,:));
                   
            % find time data point & block last byte
            temp = sum((diff(raw_data(header.timestamp-1:header.timestamp+7,:)'))>0);
            [~,tp] = max(temp);
            header.timestamp = header.timestamp-2+tp;
            raw_data(header.timestamp+7,:) = 0;
            
            epochtime(cnt_scan:cnt_scan+scan_num-1) = typecast(reshape(uint8(raw_data(header.timestamp:header.timestamp+7,:)),1,[]),'uint64');
%             header.timestamp = 16;
%             datetime(typecast(reshape(uint8(header_data(header.timestamp:header.timestamp+7,:)),1,[]),'uint64'),'ConvertFrom','epochtime','TicksPerSecond',1e3);
            raw_data(1:header.header_num,:) = [];
%             datetime(epochtime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3)
            end
            
            % raw_data
            if 1
                raw_data = raw_data(:);
                switch opt.data_ver
                    case 1
                        raw_data = reshape(raw_data,9,[]);
                        Data.EIT.OV_flag(:,cnt_scan:cnt_scan+scan_num-1) = reshape(logical(raw_data(1,:)),header.meas_num,[]);
                        Data.EIT.phase(:,cnt_scan:cnt_scan+scan_num-1) = double(reshape(typecast(reshape(uint8(raw_data(2:5,:)),1,[]),'single'),header.meas_num,[]));
                        Data.EIT.ohm(:,cnt_scan:cnt_scan+scan_num-1) = double(reshape(typecast(reshape(uint8(raw_data(6:9,:)),1,[]),'single'),header.meas_num,[]));
                    case 0
                        Data.EIT.OV_flag(:,cnt_scan:cnt_scan+scan_num-1) = reshape(raw_data(2:6:end),header.meas_num,length(raw_data)/header.meas_num/6);

                        V_R = raw_data(3:6:end) + raw_data(4:6:end).*256;
                        V_R(V_R>32767) = V_R(V_R>32767) - 65536;
                        V_I = raw_data(5:6:end) + raw_data(6:6:end).*256;
                        V_I(V_I>32767) = V_I(V_I>32767) - 65536;
                        raw_data = V_R + V_I*1i;
                        raw_data = reshape(raw_data,header.meas_num,length(raw_data)/header.meas_num);
                        switch opt.type
                            case 'magnitude'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data);
                            case 'complex'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = raw_data;
                            case 'phase'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data);
                                Data.EIT.phase(:,cnt_scan:cnt_scan+scan_num-1) = angle(raw_data);
                                if opt.mask == 1
                                    Data.EIT.phase(Data.ch_mask,:) = [];
                                end
                        end
                    case 2
                        raw_data = reshape(raw_data,256*10,[]);
                        Data.EIT.gain = reshape(typecast(uint8(reshape(raw_data(1536+1:end,:),1,[])),'single'),256,[]);
                        
                        raw_data(1536+1:end,:) = [];
                        raw_data = raw_data(:);
                        
                        Data.EIT.OV_flag(:,cnt_scan:cnt_scan+scan_num-1) = reshape(floor(raw_data(1:6:end)/(2^7)),header.meas_num,length(raw_data)/header.meas_num/6);

                        V_R = raw_data(3:6:end) + raw_data(4:6:end).*256;
                        V_R(V_R>32767) = V_R(V_R>32767) - 65536;
                        V_I = raw_data(5:6:end) + raw_data(6:6:end).*256;
                        V_I(V_I>32767) = V_I(V_I>32767) - 65536;
                        raw_data = V_R + V_I*1i;
                        raw_data = reshape(raw_data,header.meas_num,length(raw_data)/header.meas_num);
                        
                        switch opt.type
                            case 'magnitude'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data);
                            case 'complex'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = raw_data;
                            case 'phase'
                                Data.EIT.raw(:,cnt_scan:cnt_scan+scan_num-1) = abs(raw_data);
                                Data.EIT.phase(:,cnt_scan:cnt_scan+scan_num-1) = angle(raw_data);
                                if opt.mask == 1
                                    Data.EIT.phase(Data.ch_mask,:) = [];
                                end
                        end
                end
                clear raw_data;
            end
            
            % ecg
            if isfield(Path,'ecg')
                if cnt.ecg >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).ecg);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    
                    if cnt_EIT == 1
                        header.ecg_size = round(length(raw_data)/length(epochtime));
%                         header.ECG_fs = (header.ecg_size-1-4)/4; % 1:num, 4:idx
                    end
                
                    if mod(length(raw_data),37) == 0
                        header.ecg_size = 37;
                        raw_data = reshape(raw_data,header.ecg_size,[]);
                        Data.ECG.fs_ecg = header.fs_EIT;
                        temp_idx_ecg(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.ecg_idx:header.ecg_idx+3,:)),1,[]),'uint32'));
                        if median(diff(temp_idx_ecg)) == 1
                            Data.ECG.idx_ecg(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.ecg_idx:header.ecg_idx+3,:)),1,[]),'uint32'));
                            Data.ECG.ECG_hv(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = sum(raw_data(header.ecg_data:header.ecg_data+1,:).*[1 256]');
                        else
                            Data.ECG.idx_ecgx5(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.ecg_idx:header.ecg_idx+3,:)),1,[]),'uint32'));
                            Data.ECG.ECG_hvx5(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = sum(raw_data(header.ecg_data:header.ecg_data+1,:).*[1 256]');
                            Data.ECG.idx_ecg = Data.ECG.idx_ecgx5(1:5:end);
                            Data.ECG.ECG_hv = Data.ECG.ECG_hvx5(1:5:end);
                        end
                    else
                        if mod(length(raw_data),header.ecg_size) ~= 0
                            if cnt_EIT < total_EIT
                                disp('ECG data miss');
                            end
                            raw_data(end-mod(length(raw_data),header.ecg_size)+1:end) = [];
                        end
                        raw_data = reshape(raw_data,header.ecg_size,[]);
                        Data.ECG.idx_ecg(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.ecg_idx:header.ecg_idx+3,:)),1,[]),'uint32'));
                        Data.ECG.fs_ecg = header.fs_EIT;
                        Data.ECG.ECG_hv(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = sum(raw_data(header.ecg_data:header.ecg_data+1,:).*[1 256]');

                        temp = raw_data(header.ecg_data:header.ecg_data+9,:);
                        temp = sum(reshape(temp(:),2,[]).*[1 256]');
                        Data.ECG.ECG_hvx5(:,(cnt_scan-1)*5+1:(cnt_scan-1)*5+1+size(raw_data,2)*5-1) = temp;
    %                     Data.ECG.ECG_hv(:,(cnt_scan-1)*Data.ECG.ECG_fs/header.fs_EIT+1:(cnt_scan+size(raw_data,2)-1)*Data.ECG.ECG_fs/header.fs_EIT) = ...
    %                         double(typecast(reshape(uint8(raw_data(header.ecg_data:end,:)),1,[]),'uint16'));
                    end
                    clear raw_data;
                end
            end

            % Pleth
            if isfield(Path,'pleth')
                if cnt.pleth >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).pleth);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    raw_data = reshape(raw_data,9,[]);
                    
                    if ~isfield(Data,'Pleth')
                        Data.Pleth.idx_pleth = double(typecast(reshape(uint8(raw_data(2:5,:)),1,[]),'uint32'));
                        Data.Pleth.idx_Array = double(typecast(reshape(uint8(raw_data(6:9,:)),1,[]),'uint32'));
                        Data.Pleth.wave = raw_data(1,:);
                        Data.Pleth.fs_pleth = 75;
                    else
                        Data.Pleth.idx_pleth = [Data.Pleth.idx_pleth double(typecast(reshape(uint8(raw_data(2:5,:)),1,[]),'uint32'))];
                        Data.Pleth.idx_Array = [Data.Pleth.idx_Array double(typecast(reshape(uint8(raw_data(6:9,:)),1,[]),'uint32'))];
                        Data.Pleth.wave = [Data.Pleth.wave raw_data(1,:)];
                    end
                    clear raw_data;
                end
            end
            
            % SpO2
            if isfield(Path,'spO2')
                if cnt.spO2 >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).spO2);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    raw_data = reshape(raw_data,9,[]);
                    
                    if ~isfield(Data,'SpO2')
                        Data.SpO2.idx_spO2 = double(typecast(reshape(uint8(raw_data(2:5,:)),1,[]),'uint32'));
                        Data.SpO2.SpO2 = raw_data(1,:);
                        Data.SpO2.fs_spO2 = 1.5;
                    else
                        Data.SpO2.idx_spO2 = [Data.Pleth.idx_spO2 double(typecast(reshape(uint8(raw_data(2:5,:)),1,[]),'uint32'))];
                        Data.SpO2.SpO2 = [Data.Pleth.SpO2 raw_data(1,:)];
                    end
                    clear raw_data;
                end
            end
            
            % pressure
            if isfield(Path,'pressure')
                if cnt.pressure >= cnt_EIT
                    fid = fopen(Path(cnt_EIT).pressure);
                    raw_data = fread(fid,'uint8');
                    fclose(fid);
                    
                    if mod(length(raw_data),header.pres_size) ~= 0
                        if cnt_EIT < total_EIT
                            disp('Pressure data miss');
                        end
                        raw_data(end-mod(length(raw_data),header.pres_size)+1:end) = [];
                    end
                    raw_data = reshape(raw_data,header.pres_size,[]);
                    
                    header.pres_idxscan = 13:16;
                    Data.Pressure.idx_press(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = double(typecast(reshape(uint8(raw_data(header.pres_idxscan,:)),1,[]),'uint32'));
                    
                    temp_flow = raw_data(3,:) + raw_data(4,:).*256;
                    temp_flow(temp_flow>32767) = temp_flow(temp_flow>32767) - 65536;
                    Data.Pressure.Paw(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = temp_flow;
                    
                    temp_flow = raw_data(9,:) + raw_data(10,:).*256;
                    temp_flow(temp_flow>32767) = temp_flow(temp_flow>32767) - 65536;
                    Data.Pressure.Faw(:,cnt_scan:cnt_scan+size(raw_data,2)-1) = temp_flow;
                    clear raw_data;
                end
            end
            
            cnt_scan = cnt_scan + scan_num;
        end
        
        % idx checking (start number is 1 for matlab index)
        if Data.EIT.idx_scan(1) == 0
            Data.EIT.idx_scan = Data.EIT.idx_scan + 1;
        end
        
%         if sum(diff(Data.EIT.idx_scan) ~= 1)
%             if  find(diff(Data.EIT.idx_scan) ~= 1) == 1 || Data.EIT.idx_scan(2) == 1
%                 Data.EIT.idx_scan(1) = Data.EIT.idx_scan(2)-1;
%             else
%                 msgbox('scan index missing');
%             end
%         end
        
        if isfield(Path,'ecg')
            if Data.ECG.idx_ecg(1) == 0
                Data.ECG.idx_ecg = Data.ECG.idx_ecg + 1;
            end
            
            if sum(diff(Data.ECG.idx_ecg) ~= 1)
%                 if  find(diff(Data.ECG.idx_ecg) ~= 1) == 1 || Data.ECG.idx_ecg(2) == 1
%                     Data.ECG.idx_ecg(1) = Data.ECG.idx_ecg(2)-1;
%                 else
%                     msgbox('scan index missing');
%                 end
            end
        end
        
        if isfield(Path,'pressure')
            if Data.Pressure.idx_press(1) == 0
                Data.Pressure.idx_press = Data.Pressure.idx_press + 1;
            end
            
%             if sum(diff(Data.Pressure.idx_press) ~= 1)
%                 if  find(diff(Data.Pressure.idx_press) ~= 1) == 1 || Data.Pressure.idx_press(2) == 1
%                     Data.Pressure.idx_press(1) = Data.Pressure.idx_press(2)-1;
%                 else
%                     msgbox('scan index missing');
%                 end
%             end
        end
        
        % Pleth index allign
        if isfield(Path,'pleth')            
            [~,IA,~] = unique(Data.Pleth.idx_pleth,'last');
            Data.Pleth.idx_pleth_new = nan(size(Data.Pleth.idx_pleth));
            Data.Pleth.idx_pleth_new(IA) = Data.Pleth.idx_pleth(IA);
            Data.Pleth.idx_pleth_new = round(fillmissing(Data.Pleth.idx_pleth_new,'linear'));
            Data.Pleth.idx_Array(Data.Pleth.idx_pleth_new<1) = [];
            Data.Pleth.wave(Data.Pleth.idx_pleth_new<1) = [];
            Data.Pleth.idx_pleth_new(Data.Pleth.idx_pleth_new<1) = [];
            Data.Pleth.idx_pleth = Data.Pleth.idx_pleth_new;
            Data.Pleth = rmfield(Data.Pleth,'idx_pleth_new');
        end
        
        % setting file
        if isfield(Path,'setting')
            fid = fopen(Path(1).setting);
            raw_data = fread(fid,'uint8');
            fclose(fid);
            
            if rem(size(raw_data,1),(header.set_total_num)) ~= 0
                if rem(size(raw_data,1),(header.set_total_num)+8) ~= 0
                    msgbox('data missing: Setting_bin');
                    header.set_loadfail = 1;
                else
                    header.set_total_num = header.set_total_num + 8;
                    header.set_nAS = size(raw_data,1)/header.set_total_num;
                    header.set_loadfail = 0;
                    header.set_eleconfig = 1;
                end
            else
                header.set_nAS = size(raw_data,1)/header.set_total_num;
                header.set_loadfail = 0;
                header.set_eleconfig = 0;
            end
            
            if header.set_loadfail ~= 1
                raw_data = reshape(raw_data,header.set_total_num,header.set_nAS);
                
                if header.set_eleconfig == 1
                    Data.EIT.Settings.Elec_config = raw_data(1,1);
                    raw_data([1:4 end-3:end],:) = [];
                end
                
                Data.EIT.Settings.Amp = typecast(reshape(uint8(raw_data(header.set_Amp,:)),1,[]),'double');
                Data.EIT.Settings.table_CodeKey = raw_data(header.set_Gain_table_CodeKey,:);
                Data.EIT.Settings.idx_scan = raw_data(header.set_ScanIdx(4),:)*256^3 + raw_data(header.set_ScanIdx(3),:)*256^2 + raw_data(header.set_ScanIdx(2),:)*256^1 + raw_data(header.set_ScanIdx(1),:);
                
                try
                    Data.EIT.Settings.t_hms = datetime(typecast(reshape(uint8(raw_data(header.set_TimeStamp+1,:)),1,[]),'uint64'),'ConvertFrom','epochtime','TicksPerSecond',1e3);
                catch
                    try % error 'TicksPerSecond' option (matlab version)
                        Data.EIT.Settings.t_hms = datetime(typecast(reshape(uint8(raw_data(header.set_TimeStamp+1,:)),1,[]),'uint64')/1e3,'ConvertFrom','epochtime');
                    catch
                        msgbox('Datetime data is corrupted');
                    end
                end
                
                raw_data = raw_data([10:10+header.set_size_interch-1],:);
                for cnt_AS = 1:header.set_nAS
                    gain_data = reshape(raw_data(:,cnt_AS),[],16);
                    Data.EIT.Settings.ProtocolType(:,cnt_AS) = gain_data(header.set_Protocoltype,:);
                    Data.EIT.Settings.ProjNum(:,cnt_AS) = gain_data(header.set_ProjNum,:);
                    Data.EIT.Settings.Update(:,cnt_AS) = gain_data(header.set_Update,:);
                    
                    Data.EIT.Settings.GainX(:,cnt_AS) = typecast(uint8(reshape(gain_data(header.set_GainX,:),[],1)),'double');
                    Data.EIT.Settings.GainDigi1(:,cnt_AS) = reshape(gain_data(header.set_GainDigi,:),[],1);
                    Data.EIT.Settings.GainDigi2(:,cnt_AS) = reshape(gain_data(header.set_GainDigi2,:),[],1);
                end
                
                % EIT raw -> ohm
                if opt.data_ver ~= 1
                    Data.EIT.raw(:,~logical(Data.EIT.Valid)) = NaN;
                    scaler_ohm = zeros(size(Data.EIT.raw));
                    for cnt_AS = 1:header.set_nAS
%                         tp_s = Data.EIT.Settings.idx_scan(cnt_AS)+1;
                        tp_s = 1;
                        if cnt_AS ~= header.set_nAS % check next AS exist!
                            tp_e = Data.EIT.Settings.idx_scan(cnt_AS+1);
                        else % apply whole data
                            tp_e = size(Data.EIT.raw,2);
                        end
                        scaler_ohm(:,tp_s:tp_e) = repmat(header.set_scale_volt./Data.EIT.Settings.GainX(:,cnt_AS)/Data.EIT.Settings.Amp(cnt_AS),1,tp_e-tp_s+1);
                    end
                    Data.EIT.ohm = Data.EIT.raw .* scaler_ohm;
                end
            end
        end
        
        % check data masking
        if opt.data_ver == 1 && opt.mask == 1
            Data.EIT.OV_flag(Data.ch_mask,:) = [];
            Data.EIT.phase(Data.ch_mask,:) = [];
            Data.EIT.ohm(Data.ch_mask,:) = [];
        elseif opt.data_ver == 2 && opt.mask == 1
            Data.EIT.raw(Data.ch_mask,:) = [];
            Data.EIT.gain(Data.ch_mask,:) = [];
        elseif opt.data_ver == 0 && opt.mask == 1
            Data.EIT.raw(Data.ch_mask,:) = [];
            Data.EIT.ohm(Data.ch_mask,:) = [];
        end
        
        % Time convert
        try
            Data.EIT.t_hms = datetime(epochtime,'ConvertFrom','epochtime','TicksPerSecond',1e3);
        catch
            try % error 'TicksPerSecond' option (matlab version)
                Data.EIT.t_hms = datetime(epochtime/1e3,'ConvertFrom','epochtime');
            catch
                msgbox('Datetime data is corrupted');
            end
        end
        
        % ecg unit convert (->mV)
        if isfield(Path,'ecg')
            if mean(Data.ECG.ECG_hv) > 10000
                Data.ECG.ECG_hv = Data.ECG.ECG_hv/65537;
            end
            Data.ECG.ECG_hv = ((Data.ECG.ECG_hv)-15000)/3000; % unit calibration to mV
            if isfield(Data.ECG,'ECG_hvx5')
                if mean(Data.ECG.ECG_hv) > 10000
                    Data.ECG.ECG_hvx5 = Data.ECG.ECG_hvx5/65537;
                end
                Data.ECG.ECG_hvx5 = ((Data.ECG.ECG_hvx5)-15000)/3000; % unit calibration to mV
            end
        end
        
        % pressure calibration
        if isfield(Path,'pressure')
            if 0
                Data.Pressure.Paw = Data.Pressure.Paw*0.26702779+0.13351485; % Gage->Pa
                Data.Pressure.Paw = Data.Pressure.Paw*0.0101972; % Pa->cmH2O

                idx_inv = Data.Pressure.Faw>0;
                Data.Pressure.Faw = (Data.Pressure.Faw*0.00009719); % Pa->cmH2O
                Data.Pressure.Faw = -0.941*(pow2(Data.Pressure.Faw))+(52.38*abs(Data.Pressure.Faw))+1; % cmH2O->L/m
                Data.Pressure.Faw(idx_inv) = -Data.Pressure.Faw(idx_inv);
            end
            
            if 1 % new             
                Data.Pressure.Paw = 0.00267*Data.Pressure.Paw+0.0047;
                table.positive = [0,0;1000,-11.1000000000000;2000,-20.1000000000000;3000,-29.1000000000000;4000,-37.6000000000000;5000,-46.1000000000000;6000,-54.1000000000000;7000,-62.1000000000000;8000,-69.6000000000000;9000,-77.1000000000000;10000,-84.6000000000000;11000,-92.1000000000000;12000,-99.1000000000000;13000,-106.100000000000;14000,-113.100000000000;15000,-120.100000000000;16000,-127.100000000000;17000,-134.100000000000;18000,-141.100000000000;19000,-148.100000000000;20000,-155.100000000000;21000,-162.100000000000;22000,-169.100000000000;23000,-176.100000000000;24000,-183.100000000000;25000,-190.100000000000;26000,-197.100000000000;27000,-204.100000000000;28000,-211.100000000000;29000,-218.100000000000;30000,-225.100000000000;31000,-232.100000000000;32000,-239.100000000000;33000,-246.100000000000];
                table.negative = [0,0;-1000,10.5000000000000;-2000,19;-3000,27.5000000000000;-4000,36;-5000,44.1000000000000;-6000,51.6000000000000;-7000,59.1000000000000;-8000,65.6000000000000;-9000,72.1000000000000;-10000,78.1000000000000;-11000,84.1000000000000;-12000,90.1000000000000;-13000,96.1000000000000;-14000,102.100000000000;-15000,108.100000000000;-16000,114.100000000000;-17000,120.100000000000;-18000,126.100000000000;-19000,132.100000000000;-20000,138.100000000000;-21000,144.100000000000;-22000,150.100000000000;-23000,156.100000000000;-24000,162.100000000000;-25000,168.100000000000;-26000,174.100000000000;-27000,180.100000000000;-28000,186.100000000000;-29000,192.100000000000;-30000,198.100000000000;-31000,204.100000000000;-32000,210.100000000000;-33000,216.100000000000];
                temp = Data.Pressure.Faw;
                idx_positive=Data.Pressure.Faw>0;
                temp_p = temp(idx_positive);
                temp(idx_positive) = interp1(table.positive(:,1),table.positive(:,2),temp_p);
                temp_n = temp(~idx_positive);
                temp(~idx_positive) = interp1(table.negative(:,1),table.negative(:,2),temp_n);
                Data.Pressure.Faw = temp;
            end
            
            if Data.Pressure.idx_press(1) == 1 % first data error
                Data.Pressure.idx_press = [0 Data.Pressure.idx_press];
                Data.Pressure.Paw = [0 Data.Pressure.Paw];
                Data.Pressure.Faw = [0 Data.Pressure.Faw];
            end
            
%             if length(Data.Pressure.idx_press) ~= length(epochtime) % data number missmatch
%                 if length(Data.Pressure.idx_press) < length(epochtime)
%                     Data.Pressure.idx_press(length(Data.Pressure.idx_press)+1:length(epochtime)) = 0;
%                     Data.Pressure.Paw(length(Data.Pressure.idx_press)+1:length(epochtime)) = 0;
%                     Data.Pressure.Faw(length(Data.Pressure.idx_press)+1:length(epochtime)) = 0;
%                 else
%                     Data.Pressure.idx_press(length(epochtime)+1:end) = [];
%                     Data.Pressure.Paw(length(epochtime)+1:end) = [];
%                     Data.Pressure.Faw(length(epochtime)+1:end) = [];
%                 end
%             end
        end
        
%         if sum(diff(epochtime)<=0) > 1
%             rm_idx = [false diff(epochtime)<=0];
%             FN = fieldnames(Data.EIT);
%             for cnt_field = 1:length(FN)
%                 if ~isstruct(Data.EIT.(FN{cnt_field}))
%                     Data.EIT.(FN{cnt_field})(:,rm_idx) = [];
%                 end
%             end
%             if isfield(Path,'pressure')
%                 Data.Pressure.idx_press(rm_idx) = [];
% %                 FN = fieldnames(Data.Pressure);
% %                 for cnt_field = 1:length(FN)
% %                     if ~isstruct(Data.Pressure.(FN{cnt_field}))
% %                         Data.Pressure.(FN{cnt_field})(rm_idx) = [];
% %                     end
% %                 end
%             end
%         end 
end
    if opt.disp
        delete(f);
    end
catch ME    
    msgbox('EIT Import failed');
    if opt.disp
        delete(f);
    end
    rethrow(ME)
end
    
end

function [mask] = FxEIT_mask(ch,first_sat)
    temp = zeros(ch,ch);
    if nargin == 1
        first_sat = [ch,1,2];
    end

    if length(first_sat) == 4
        cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3); cnt4 = first_sat(4);
        for i = 1:ch
            temp([cnt1,cnt2,cnt3,cnt4],i) = 1;
            cnt1 = cnt1 + 1;
            cnt2 = cnt2 + 1;
            cnt3 = cnt3 + 1;
            cnt4 = cnt4 + 1;
            if cnt1 > ch
                cnt1 = 1;
            end
            if cnt2 > ch
                cnt2 = 1;
            end
            if cnt3 > ch
                cnt3 = 1;
            end
            if cnt4 > ch
                cnt4 = 1;
            end
        end
    else
        cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3);
        for i = 1:ch
            temp([cnt1,cnt2,cnt3],i) = 1;
            cnt1 = cnt1 + 1;
            cnt2 = cnt2 + 1;
            cnt3 = cnt3 + 1;
            if cnt1 > ch
                cnt1 = 1;
            end
            if cnt2 > ch
                cnt2 = 1;
            end
            if cnt3 > ch
                cnt3 = 1;
            end
        end
    end

    temp2 = reshape(temp,ch*ch,1);
    [mask, ~] = find(temp2 == 1);
       
end

 
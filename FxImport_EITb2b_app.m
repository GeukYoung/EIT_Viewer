function [Data] = FxImport_EITb2b_app(path)
if nargin < 1
    path{1} = uigetdir; 
end
Data = [];
% path{1} = 'D:\OneDrive\2. Project\01. Research\Bilab work\7. newTVData\Test_20220825_161018';
% path{1} = 'D:\OneDrive\2. Project\01. Research\Bilab work\7. newTVData\LongTerm_20220825_170121\Parameter';
% path{1} = 'D:\OneDrive\2. Project\01. Research\Bilab work\7. newTVData\20220829_addROI';

% 폴더형식 미정
if ischar(path)
    Path_main = path;
    if ~contains(Path_main, '\Parameter')
        temp_dir = dir(Path_main);
        if contains([temp_dir.name], 'Parameter')
            Path_main = [Path_main '\Parameter'];
        end
        clear temp_dir;
    end
%     path = path(2:end);
%     argOffset = 1;
else
    error('ERROR: Data path');
end

%% call file list
dirlist = dir(Path_main);
cnt.Wave = 0; % AT & HV
cnt.ScaleFactor = 0; % AT & HV

cnt.PFV = 0; % AT
cnt.EIT_AT = 0; % AT
cnt.ROI = 0; % AT

cnt.CVSTime = 0; % HT
cnt.RWaveTime = 0; % HT
cnt.SpO2 = 0; % HT
cnt.Z = 0; % HT
cnt.NIBP = 0; % HT
cnt.SV = 0; % HT
file_Type = '';

% check file contents
for cnt_file = 1:length(dirlist)
    temp_string = strsplit(dirlist(cnt_file).name,{'.', '['});
    idx.file = str2double(regexp(temp_string{1},'\d*','match'));
    if ~isempty(idx.file)
        idx.file = idx.file(end) + 1;
    end

    if contains(temp_string{1},'GraphAndSensor_','IgnoreCase',true)
        idx.Wave(idx.file) = true;
        cnt.Wave = cnt.Wave + 1;
        Path(idx.file).Wave = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'ScaleFactor_','IgnoreCase',true)
        idx.ScaleFactor(idx.file) = true;
        cnt.ScaleFactor = cnt.ScaleFactor + 1;
        Path(idx.file).ScaleFactor = fullfile(Path_main, dirlist(cnt_file).name);
    end
 
    if contains(temp_string{1},'RVDEIT_','IgnoreCase',true)
        idx.EIT(idx.file) = true;
        cnt.EIT_AT = cnt.EIT_AT + 1;
        Path(idx.file).EIT = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'PFV_','IgnoreCase',true)
        idx.PFV(idx.file) = true;
        cnt.PFV = cnt.PFV + 1;
        Path(idx.file).PFV = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'ROI_','IgnoreCase',true)
        idx.ROI(idx.file) = true;
        cnt.ROI = cnt.ROI + 1;
        Path(idx.file).ROI = fullfile(Path_main, dirlist(cnt_file).name);
    end

    if contains(temp_string{1},'CVSTime_','IgnoreCase',true)
        idx.CVSTime(idx.file) = true;
        cnt.CVSTime = cnt.CVSTime + 1;
        Path(idx.file).CVSTime = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'RWaveTime_','IgnoreCase',true)
        idx.RWaveTime(idx.file) = true;
        cnt.RWaveTime = cnt.RWaveTime + 1;
        Path(idx.file).RWaveTime = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'SpO2Relate_','IgnoreCase',true)
        idx.SpO2(idx.file) = true;
        cnt.SpO2 = cnt.SpO2 + 1;
        Path(idx.file).SpO2 = fullfile(Path_main, dirlist(cnt_file).name);
    end
    
    if contains(temp_string{1},'Z_','IgnoreCase',true)
        idx.Z(idx.file) = true;
        cnt.Z = cnt.Z + 1;
        Path(idx.file).Z = fullfile(Path_main, dirlist(cnt_file).name);    
    end
    
    if contains(temp_string{1},'NIBP_','IgnoreCase',true)
        idx.NIBP(idx.file) = true;
        cnt.NIBP = cnt.NIBP + 1;
        Path(idx.file).NIBP = fullfile(Path_main, dirlist(cnt_file).name);    
    end
    
    if contains(temp_string{1},'SV_','IgnoreCase',true)
        idx.SV(idx.file) = true;
        cnt.SV = cnt.SV + 1;
        Path(idx.file).SV = fullfile(Path_main, dirlist(cnt_file).name);    
    end
end

if cnt.EIT_AT > 0
    file_Type = 'AT';
else
    file_Type = 'HV';
end

%% ScaleFactor
if cnt.ScaleFactor~=0
    for cnt_file = 1:cnt.ScaleFactor
       [temp{cnt_file}, header] = Import_ScaleFactor(Path(cnt_file).ScaleFactor);
    end
    
    % add file_Type to separate (TV/SV scalefactor)
    %
    % 
    %
    
    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.ScaleFactor = structcat(temp);
    clear temp;
end

%% Waveform (GraphAndSensor)
if cnt.Wave~=0
    for cnt_file = 1:cnt.Wave
        switch file_Type
            case 'AT'
                [temp{cnt_file}, header] = Import_Wave_AT(Path(cnt_file).Wave);
            case 'HV'
                [temp{cnt_file}, header] = Import_Wave_HV(Path(cnt_file).Wave);
        end
    end
    
    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.Wave = structcat(temp);
    clear temp;
end

%% PFV
if cnt.PFV~=0
    for cnt_file = 1:cnt.PFV
        [temp{cnt_file}, header] = Import_PFV(Path(cnt_file).PFV);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.PFV = structcat(temp);
    clear temp;
end

%% AT param (RVDEIT)
if cnt.EIT_AT~=0
    for cnt_file = 1:cnt.EIT_AT
        [temp{cnt_file}, header] = Import_EIT(Path(cnt_file).EIT);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.EIT = structcat(temp);
    clear temp;
end

%% CVS peak (CVSTime)
if cnt.CVSTime~=0
    for cnt_file = 1:cnt.CVSTime
        [temp{cnt_file}, header] = Import_CVS_peakvalley(Path(cnt_file).CVSTime);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.peak_CVS = structcat(temp);
    clear temp;
end
 
%% ECG peak (RWaveTime)
if cnt.RWaveTime~=0
    for cnt_file = 1:cnt.RWaveTime
        [temp{cnt_file}, header] = Import_ECG_peakvalley(Path(cnt_file).RWaveTime);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.peak_ECG = structcat(temp);
    clear temp;
end

%% SpO2
if cnt.SpO2~=0
    for cnt_file = 1:cnt.SpO2
        [temp{cnt_file}, header] = Import_SpO2(Path(cnt_file).SpO2);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.SpO2 = structcat(temp);
    clear temp;
end

%% Z_0
if cnt.Z~=0
    for cnt_file = 1:cnt.Z
        [temp{cnt_file}, header] = Import_Z_0(Path(cnt_file).Z);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.Z = structcat(temp);
    clear temp;
end

%% NIBP
if cnt.NIBP~=0
    for cnt_file = 1:cnt.NIBP
        [temp{cnt_file}, header] = Import_NIBP(Path(cnt_file).NIBP);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.NIBP = structcat(temp);
    clear temp;
end

%% SV
if cnt.SV~=0
    for cnt_file = 1:cnt.SV
        [temp{cnt_file}, header] = Import_SV(Path(cnt_file).SV);
    end

    if isfield(Data,'Header') ~= 1
        Data.Header = header;
    else
        Data.Header = appendStruct(Data.Header, header);
    end
    Data.SV = structcat(temp);
    clear temp;
end
end

function [Data, Header] = Import_Wave_AT(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_RVS = 53:56;
    instruction.ParameterDataType_SQI_RVS = 57:60;
    instruction.ParameterDataType_Paw = 61:64;
    instruction.ParameterDataType_Faw = 65:68;
    instruction.ParameterDataType_Vaw = 69:72;
    instruction.patientID_Length = 73:76;
    instruction.patientID_end = double(77+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 77:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.RVS = 13:16;
    instruction.SQI = 17:20;
    instruction.Paw = 21:24;
    instruction.Faw = 25:28;
    instruction.Vaw = 29:32;
    instruction.nData = 32;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_RVS = typecast(raw_data(instruction.ParameterDataType_RVS,:),'uint32');
    Header.ParameterDataType_SQI_RVS = typecast(raw_data(instruction.ParameterDataType_SQI_RVS,:),'uint32');
    Header.ParameterDataType_Paw = typecast(raw_data(instruction.ParameterDataType_Paw,:),'uint32');
    Header.ParameterDataType_Faw = typecast(raw_data(instruction.ParameterDataType_Faw,:),'uint32');
    Header.ParameterDataType_Vaw = typecast(raw_data(instruction.ParameterDataType_Vaw,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.RVS = typecast(uint8(reshape(raw_data(instruction.RVS,:),[],1)),'single')';
    Data.SQI = typecast(uint8(reshape(raw_data(instruction.SQI,:),[],1)),'uint32')';
    Data.Paw = typecast(uint8(reshape(raw_data(instruction.Paw,:),[],1)),'single')';
    Data.Faw = typecast(uint8(reshape(raw_data(instruction.Faw,:),[],1)),'single')';
    Data.Vaw = typecast(uint8(reshape(raw_data(instruction.Vaw,:),[],1)),'single')';
end

function [Data, Header] = Import_Wave_HV(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_CVS = 53:56;
    instruction.ParameterDataType_SQI_CVS = 57:60;
    instruction.ParameterDataType_ECG = 61:64;
    instruction.patientID_Length = 65:68;
    instruction.patientID_end = double(69+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 69:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.CVS = 13:16;
    instruction.SQI = 17:20;
    instruction.ECG = 21:30;
    instruction.nData = 30;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_CVS = typecast(raw_data(instruction.ParameterDataType_CVS,:),'uint32');
    Header.ParameterDataType_SQI_CVS = typecast(raw_data(instruction.ParameterDataType_SQI_CVS,:),'uint32');
    Header.ParameterDataType_ECG = typecast(raw_data(instruction.ParameterDataType_ECG,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,30,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.CVS = typecast(uint8(reshape(raw_data(instruction.CVS,:),[],1)),'single')';
    Data.SQI = typecast(uint8(reshape(raw_data(instruction.SQI,:),[],1)),'uint32')';
    Data.ECG = typecast(uint8(reshape(raw_data(instruction.ECG,:),[],1)),'uint16')';
    Data.ECG = [Data.ECG(46:50:end); Data.ECG(47:50:end); Data.ECG(48:50:end); Data.ECG(49:50:end); Data.ECG(50:50:end)];
end

function [Data, Header] = Import_PFV(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_PEEP = 53:56;
    instruction.ParameterDataType_tPEEP = 57:60;
    instruction.ParameterDataType_Ppeak = 61:64;
    instruction.ParameterDataType_tPpeak = 65:68;
    instruction.ParameterDataType_Pplat = 69:72;
    instruction.ParameterDataType_tPplat = 73:76;
    instruction.ParameterDataType_Fmax = 77:80;
    instruction.ParameterDataType_tFmax = 81:84;
    instruction.ParameterDataType_Fmin = 85:88;
    instruction.ParameterDataType_tFmin = 89:92;
    instruction.ParameterDataType_TVi = 93:96;
    instruction.ParameterDataType_tTVi = 97:100;
    instruction.ParameterDataType_TVe = 101:104;
    instruction.ParameterDataType_tTVe = 105:108;
    instruction.patientID_Length = 109:112;
    instruction.patientID_end = double(113+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 113:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.PEEP = 13:16;
    instruction.tPEEP = 17:24;
    instruction.Ppeak = 25:28;
    instruction.tPpeak = 29:36;
    instruction.Pplat = 37:40;
    instruction.tPplat = 41:48;
    instruction.Fmax = 49:52;
    instruction.tFmax = 53:60;
    instruction.Fmin = 61:64;
    instruction.tFmin = 65:72;
    instruction.TVi = 73:76;
    instruction.tTVi = 77:84;
    instruction.TVe = 85:88;
    instruction.tTVe = 89:96;
    instruction.nData = 96;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_PEEP = typecast(raw_data(instruction.ParameterDataType_PEEP,:),'uint32');
    Header.ParameterDataType_tPEEP = typecast(raw_data(instruction.ParameterDataType_tPEEP,:),'uint32');
    Header.ParameterDataType_Ppeak = typecast(raw_data(instruction.ParameterDataType_Ppeak,:),'uint32');
    Header.ParameterDataType_tPpeak = typecast(raw_data(instruction.ParameterDataType_tPpeak,:),'uint32');
    Header.ParameterDataType_Pplat = typecast(raw_data(instruction.ParameterDataType_Pplat,:),'uint32');
    Header.ParameterDataType_tPplat = typecast(raw_data(instruction.ParameterDataType_tPplat,:),'uint32');
    Header.ParameterDataType_Fmax = typecast(raw_data(instruction.ParameterDataType_Fmax,:),'uint32');
    Header.ParameterDataType_tFmax = typecast(raw_data(instruction.ParameterDataType_tFmax,:),'uint32');
    Header.ParameterDataType_Fmin = typecast(raw_data(instruction.ParameterDataType_Fmin,:),'uint32');
    Header.ParameterDataType_tFmin = typecast(raw_data(instruction.ParameterDataType_tFmin,:),'uint32');
    Header.ParameterDataType_TVi = typecast(raw_data(instruction.ParameterDataType_TVi,:),'uint32');
    Header.ParameterDataType_tTVi = typecast(raw_data(instruction.ParameterDataType_tTVi,:),'uint32');
    Header.ParameterDataType_TVe = typecast(raw_data(instruction.ParameterDataType_TVe,:),'uint32');
    Header.ParameterDataType_tTVe = typecast(raw_data(instruction.ParameterDataType_tTVe,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.PEEP = typecast(uint8(reshape(raw_data(instruction.PEEP,:),[],1)),'single')';
    Data.tPEEP = typecast(uint8(reshape(raw_data(instruction.tPEEP,:),[],1)),'uint64')';
    Data.tPEEP = datetime(Data.tPEEP,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.Ppeak = typecast(uint8(reshape(raw_data(instruction.Ppeak,:),[],1)),'single')';
    Data.tPpeak = typecast(uint8(reshape(raw_data(instruction.tPpeak,:),[],1)),'uint64')';
    Data.tPpeak = datetime(Data.tPpeak,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.Pplat = typecast(uint8(reshape(raw_data(instruction.Pplat,:),[],1)),'single')';
    Data.tPplat = typecast(uint8(reshape(raw_data(instruction.tPplat,:),[],1)),'uint64')';
    Data.tPplat = datetime(Data.tPplat,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.Fmax = typecast(uint8(reshape(raw_data(instruction.Fmax,:),[],1)),'single')';
    Data.tFmax = typecast(uint8(reshape(raw_data(instruction.tFmax,:),[],1)),'uint64')';
    Data.tFmax = datetime(Data.tFmax,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.Fmin = typecast(uint8(reshape(raw_data(instruction.Fmin,:),[],1)),'single')';
    Data.tFmin = typecast(uint8(reshape(raw_data(instruction.tFmin,:),[],1)),'uint64')';
    Data.tFmin = datetime(Data.tFmin,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.TVi = typecast(uint8(reshape(raw_data(instruction.TVi,:),[],1)),'single')';
    Data.tTVi = typecast(uint8(reshape(raw_data(instruction.tTVi,:),[],1)),'uint64')';
    Data.tTVi = datetime(Data.tTVi,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.TVe = typecast(uint8(reshape(raw_data(instruction.TVe,:),[],1)),'single')';
    Data.tTVe = typecast(uint8(reshape(raw_data(instruction.tTVe,:),[],1)),'uint64')';
    Data.tTVe = datetime(Data.tTVe,'ConvertFrom','epochtime','TicksPerSecond',1e3);
end

function [Data, Header] = Import_EIT(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_RVSvalleyTime = 53:56;
    instruction.ParameterDataType_RVSpeakTime = 57:60;
    instruction.ParameterDataType_EITvalley = 61:64;
    instruction.ParameterDataType_EITpeak = 65:68;
    instruction.ParameterDataType_RVD = 69:72;
    instruction.patientID_Length = 73:76;
    instruction.patientID_end = double(77+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 77:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.tRVSvalley = 13:20;
    instruction.tRVSpeak = 21:28;
    instruction.EITvalley_Z = 29:(29+1664-1); % 8(Zi+Zr) x 208
    instruction.EITpeak_Z = 1693:(1693+1664-1);
    instruction.nEIT = 208;
    instruction.ImRVD = 3357:19740;
    instruction.nRVD = 64*64;
    instruction.nData = 19740;
    
%     instruction.EITvalley = 29:2356; % 24(header) + 9(OV+Zi+Zr) x 256
%     instruction.EITvalley_h = instruction.EITvalley(1:24); instruction.EITvalley(1:24)=[];
%     instruction.EITvalley_OV = instruction.EITvalley(1:9:end); instruction.EITvalley(1:9:end) = [];
%     instruction.EITvalley_Z = instruction.EITvalley;
%     instruction.EITpeak = 2357:4684;
%     instruction.EITpeak_h = instruction.EITpeak(1:24); instruction.EITpeak(1:24)=[];
%     instruction.EITpeak_OV = instruction.EITpeak(1:9:end); instruction.EITpeak(1:9:end) = [];
%     instruction.EITpeak_Z = instruction.EITpeak;
%     instruction.nEIT = 256;
%     instruction.ImRVD = 4685:21068;
%     instruction.nRVD = 64*64;
%     instruction.nData = 21068;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_RVSvalleyTime = typecast(raw_data(instruction.ParameterDataType_RVSvalleyTime,:),'uint32');
    Header.ParameterDataType_RVSpeakTime = typecast(raw_data(instruction.ParameterDataType_RVSpeakTime,:),'uint32');
    Header.ParameterDataType_EITvalley = typecast(raw_data(instruction.ParameterDataType_EITvalley,:),'uint32');
    Header.ParameterDataType_EITpeak = typecast(raw_data(instruction.ParameterDataType_EITpeak,:),'uint32');
    Header.ParameterDataType_RVD = typecast(raw_data(instruction.ParameterDataType_RVD,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]); %21068  10534x2
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.tRVSvalley = typecast(uint8(reshape(raw_data(instruction.tRVSvalley,:),[],1)),'uint64')';
    Data.tRVSvalley = datetime(Data.tRVSvalley,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.tRVSpeak = typecast(uint8(reshape(raw_data(instruction.tRVSpeak,:),[],1)),'uint64')';
    Data.tRVSpeak = datetime(Data.tRVSpeak,'ConvertFrom','epochtime','TicksPerSecond',1e3);

%     Data.OV_valley = reshape(typecast(uint8(reshape(raw_data(instruction.EITvalley_OV,:),[],1)),'uint8'),instruction.nEIT,[]);
    Data.EITvalley = reshape(typecast(uint8(reshape(raw_data(instruction.EITvalley_Z,:),[],1)),'single'),instruction.nEIT*2,[]);
    [Data.Zphase_valley, Data.Zmag_valley] = deal(Data.EITvalley(1:2:end,:), Data.EITvalley(2:2:end,:));
    Data = rmfield(Data,'EITvalley')';

%     Data.OV_peak = reshape(typecast(uint8(reshape(raw_data(instruction.EITpeak_OV,:),[],1)),'uint8'),instruction.nEIT,[]);
    Data.EITpeak = reshape(typecast(uint8(reshape(raw_data(instruction.EITpeak_Z,:),[],1)),'single'),instruction.nEIT*2,[]);
    [Data.Zphase_peak, Data.Zmag_peak] = deal(Data.EITpeak(1:2:end,:), Data.EITpeak(2:2:end,:));
    Data = rmfield(Data,'EITpeak')';

    Data.Im_RVD = reshape(typecast(uint8(reshape(raw_data(instruction.ImRVD,:),[],1)),'single'),instruction.nRVD,[]);
    Data.Im_RVD(Data.Im_RVD == -1000000) = NaN; % figure; imagesc(reshape(DataSet.ImRVD(:,100),64,64));
end

function [Data, Header] = Import_ScaleFactor(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_ScaleFactorTV = 53:56;
    instruction.ParameterDataType_GainTable = 57:60;
    instruction.ParameterDataType_CableDirection = 61:64;
    instruction.ParameterDataType_Amplitude = 65:68;
    instruction.patientID_Length = 69:72;
    instruction.patientID_end = double(73+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 73:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4; % int
    instruction.TimeStamp = 5:12; % uint64
    instruction.ScaleFactor = 13:16; % float
    instruction.GainTable = 17:2625; %
    instruction.CableDirection = 2626:2629; % int
    instruction.Amplitude = 2630:2633; % float
    instruction.nData = 2633;
    
    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_ScaleFactorTV = typecast(raw_data(instruction.ParameterDataType_ScaleFactorTV,:),'uint32');
    Header.ParameterDataType_GainTable = typecast(raw_data(instruction.ParameterDataType_GainTable,:),'uint32');
    Header.ParameterDataType_CableDirection = typecast(raw_data(instruction.ParameterDataType_CableDirection,:),'uint32');
    Header.ParameterDataType_Amplitude = typecast(raw_data(instruction.ParameterDataType_Amplitude,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.ScaleFactor = typecast(uint8(reshape(raw_data(instruction.ScaleFactor,:),[],1)),'single')';
%     DataSet.GainTable = typecast(uint8(reshape(raw_data(header.GainTable,:),[],1)),'uint32')';
    Data.GainTable = uint8(reshape(raw_data(instruction.GainTable,:),[],1))';
    Data.CableDirection = typecast(uint8(reshape(raw_data(instruction.CableDirection,:),[],1)),'uint32')';
    Data.Amplitude = typecast(uint8(reshape(raw_data(instruction.Amplitude,:),[],1)),'single')';
end

function [Data, Header] = Import_ROI(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_ROI = 53:56;
    instruction.patientID_Length = 57:60;
    instruction.patientID_end = double(61+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 61:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4; % int
    instruction.TimeStamp = 5:12; % uint64
    instruction.ROI0_x = 13:16; % float
    instruction.ROI0_y = 17:20; % float
    instruction.ROI0_width = 21:24; % float
    instruction.ROI0_height = 25:28; % float
    instruction.ROI1_x = 29:32; % float
    instruction.ROI1_y = 33:36; % float
    instruction.ROI1_width = 37:40; % float
    instruction.ROI1_height = 41:44; % float
    instruction.ROI2_x = 45:48; % float
    instruction.ROI2_y = 49:52; % float
    instruction.ROI2_width = 53:56; % float
    instruction.ROI2_height = 57:60; % float
    instruction.ROI3_x = 61:64; % float
    instruction.ROI3_y = 65:68; % float
    instruction.ROI3_width = 69:72; % float
    instruction.ROI3_height = 73:76; % float
    
    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_ROI = typecast(raw_data(instruction.ParameterDataType_ROI,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.ROI_1 = [typecast(uint8(reshape(raw_data(instruction.ROI0_x,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI0_y,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI0_width,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI0_height,:),[],1)),'single')'];
    Data.ROI_2 = [typecast(uint8(reshape(raw_data(instruction.ROI1_x,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI1_y,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI1_width,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI1_height,:),[],1)),'single')'];
    Data.ROI_3 = [typecast(uint8(reshape(raw_data(instruction.ROI2_x,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI2_y,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI2_width,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI2_height,:),[],1)),'single')'];
    Data.ROI_4 = [typecast(uint8(reshape(raw_data(instruction.ROI3_x,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI3_y,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI3_width,:),[],1)),'single')', ...
                  typecast(uint8(reshape(raw_data(instruction.ROI3_height,:),[],1)),'single')'];
end

function [Data, Header] = Import_CVS_peakvalley(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_CVSvalleyTime = 53:56;
    instruction.ParameterDataType_CVSpeakTime = 57:60;
    instruction.patientID_Length = 61:64;
    instruction.patientID_end = double(65+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 65:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.tValley_CVS = 13:20;
    instruction.tPeak_CVS = 21:28;
    instruction.nData = 28;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_CVSvalleyTime = typecast(raw_data(instruction.ParameterDataType_CVSvalleyTime,:),'uint32');
    Header.ParameterDataType_CVSpeakTime = typecast(raw_data(instruction.ParameterDataType_CVSpeakTime,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.tValley_CVS = typecast(uint8(reshape(raw_data(instruction.tValley_CVS,:),[],1)),'uint64')';
    Data.tValley_CVS = datetime(Data.tValley_CVS,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.tPeak_CVS = typecast(uint8(reshape(raw_data(instruction.tPeak_CVS,:),[],1)),'uint64')';
    Data.tPeak_CVS = datetime(Data.tPeak_CVS,'ConvertFrom','epochtime','TicksPerSecond',1e3);
end

function [Data, Header] = Import_ECG_peakvalley(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_RWaveTime = 53:56;
    instruction.patientID_Length = 57:60;
    instruction.patientID_end = double(61+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 61:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.tPeak_ECG = 13:20;
    instruction.nData = 20;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_RWaveTime = typecast(raw_data(instruction.ParameterDataType_RWaveTime,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.tPeak_ECG = typecast(uint8(reshape(raw_data(instruction.tPeak_ECG,:),[],1)),'uint64')';
    Data.tPeak_ECG = datetime(Data.tPeak_ECG,'ConvertFrom','epochtime','TicksPerSecond',1e3);
end

function [Data, Header] = Import_SpO2(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_Pleth = 53:56;
    instruction.ParameterDataType_SpO2 = 57:60;
    instruction.patientID_Length = 61:64;
    instruction.patientID_end = double(65+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 65:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.Pleth = 13:37;
    instruction.SpO2 = 38;
    instruction.PR_SpO2 = 39:40;
    instruction.nData = 40;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_Pleth = typecast(raw_data(instruction.ParameterDataType_Pleth,:),'uint32');
    Header.ParameterDataType_SpO2 = typecast(raw_data(instruction.ParameterDataType_SpO2,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.Pleth = uint8(reshape(raw_data(instruction.Pleth,:),[],1))';
    Data.SpO2 = int8(reshape(raw_data(instruction.SpO2,:),[],1))';
    Data.PR_SpO2 = typecast(uint8(reshape(raw_data(instruction.PR_SpO2,:),[],1)),'uint16')';
end

function [Data, Header] = Import_Z_0(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_Z = 53:56;
    instruction.patientID_Length = 57:60;
    instruction.patientID_end = double(61+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 61:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.Z = 13:16;
    instruction.nData = 16;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_Z = typecast(raw_data(instruction.ParameterDataType_Z,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.Z = typecast(uint8(reshape(raw_data(instruction.Z,:),[],1)),'single')';
end

function [Data, Header] = Import_NIBP(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_SBP = 53:56;
    instruction.ParameterDataType_DBP = 57:60;
    instruction.ParameterDataType_MAP = 61:64;
    instruction.ParameterDataType_PulseRateNIBP = 65:68;
    instruction.patientID_Length = 69:72;
    instruction.patientID_end = double(73+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 73:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.SBP = 13:14;
    instruction.DBP = 15:16;
    instruction.MAP = 17:18;
    instruction.PR_NIBP = 19:20;
    instruction.nData = 20;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_SBP = typecast(raw_data(instruction.ParameterDataType_SBP,:),'uint32');
    Header.ParameterDataType_DBP = typecast(raw_data(instruction.ParameterDataType_DBP,:),'uint32');
    Header.ParameterDataType_MAP = typecast(raw_data(instruction.ParameterDataType_MAP,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.SBP = typecast(uint8(reshape(raw_data(instruction.SBP,:),[],1)),'uint16')';
    Data.DBP = typecast(uint8(reshape(raw_data(instruction.DBP,:),[],1)),'uint16')';
    Data.MAP = typecast(uint8(reshape(raw_data(instruction.MAP,:),[],1)),'uint16')';
    Data.PR_NIBP = typecast(uint8(reshape(raw_data(instruction.PR_NIBP,:),[],1)),'uint16')';
end

function [Data, Header] = Import_SV(Path)
    % data import
    fid = fopen(Path);
    raw_data = uint8(fread(fid,'uint8'));
    fclose(fid);

    % Instruction_Header
    instruction.fileType = 1:4;
    instruction.productType = 5:8;
    instruction.headerSize = 9:12;
    instruction.headerCheckSum = 13:16;
    instruction.initialTimeStamp = 17:24;
    instruction.refTimeStamp = 25:32;
    instruction.step = 33:36;
    instruction.dataFormatVersion = 37:40;
    instruction.delimiterKey1 = 41:44;
    instruction.param_type = 45:48;
    instruction.numberOfParameter = 49:52;
    instruction.ParameterDataType_SV = 53:56;
    instruction.patientID_Length = 57:60;
    instruction.patientID_end = double(61+typecast(raw_data(instruction.patientID_Length,:),'uint32')*4-1);
    instruction.id_char = 61:instruction.patientID_end;
    instruction.patient_gender = (1:4)+instruction.patientID_end;
    instruction.patient_age = (5:8)+instruction.patientID_end;
    instruction.patient_height = (9:12)+instruction.patientID_end;
    instruction.patient_weight = (13:16)+instruction.patientID_end;

    % Instruction_Data
    instruction.delimiterKey2 = 1:4;
    instruction.TimeStamp = 5:12;
    instruction.SV = 13:16;
    instruction.nData = 16;

    % Header
    Header.fileType = char(raw_data(instruction.fileType)');
    Header.productType = typecast(raw_data(instruction.productType,:),'uint32');
    Header.headerSize = typecast(raw_data(instruction.headerSize,:),'uint32');
    Header.headerCheckSum = typecast(raw_data(instruction.headerCheckSum,:),'uint32');
    Header.initialTimeStamp = typecast(raw_data(instruction.initialTimeStamp,:),'uint64');
    Header.initialTimeStamp = datetime(Header.initialTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.refTimeStamp = typecast(raw_data(instruction.refTimeStamp,:),'uint64');
    Header.refTimeStamp = datetime(Header.refTimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Header.step = typecast(raw_data(instruction.step,:),'uint32');
    Header.dataFormatVersion = typecast(raw_data(instruction.dataFormatVersion,:),'uint32');
    Header.delimiterKey = typecast(raw_data(instruction.delimiterKey1,:),'uint32');
    Header.param_type = typecast(raw_data(instruction.param_type,:),'uint32');
    Header.numberOfParameter = typecast(raw_data(instruction.numberOfParameter,:),'uint32');
    Header.ParameterDataType_SV = typecast(raw_data(instruction.ParameterDataType_SV,:),'uint32');
    Header.patientID_Length = typecast(raw_data(instruction.patientID_Length,:),'uint32');
    Header.id_char = char(raw_data(instruction.id_char)'); % int(?) char(?)
    Header.patient_gender = typecast(raw_data(instruction.patient_gender,:),'uint32');
    Header.patient_age = typecast(raw_data(instruction.patient_age,:),'uint32');
    Header.patient_height = typecast(raw_data(instruction.patient_height,:),'single');
    Header.patient_weight = typecast(raw_data(instruction.patient_weight,:),'single');
    raw_data(1:Header.headerSize) = [];

    % Data
    raw_data = reshape(raw_data,instruction.nData,[]);
    Data.delimiterKey2 = typecast(uint8(reshape(raw_data(instruction.headerSize,:),[],1)),'uint32')';
    Data.TimeStamp = typecast(uint8(reshape(raw_data(instruction.TimeStamp,:),[],1)),'uint64')';
    Data.TimeStamp = datetime(Data.TimeStamp,'ConvertFrom','epochtime','TicksPerSecond',1e3);
    Data.SV = typecast(uint8(reshape(raw_data(instruction.SV,:),[],1)),'single')';
end

function Y = structcat(inputs)
N = length(inputs);

FN = cell(N,1) ;
VAL = cell(N,1) ;
% parse the inputs
for cnt_s = 1:N
    X = inputs{cnt_s} ;
    if ~isstruct(X)
        error('catstruct:InvalidArgument',['Argument #' num2str(cnt_s) ' is not a structure.']) ;
    end
    
    FN{cnt_s} = fieldnames(X);
    VAL{cnt_s} = struct2cell(X);
end

for cnt_f = 1:length(VAL{1})
    if ~isstruct(VAL{1}{cnt_f})
        temp = [];
        for cnt_s = 1:N
            if size(VAL{cnt_s}{cnt_f},1) == 1
                temp = [temp, VAL{cnt_s}{cnt_f}];
            else
                temp = [temp, VAL{cnt_s}{cnt_f}];
            end
        end
        Y.(FN{1}{cnt_f}) = temp;
    else
        for cnt_s = 1:N
            VAL_iner{:,cnt_s} = struct2cell(VAL{cnt_s}{cnt_f});
            FN_iner{:,cnt_s} = fieldnames(VAL{cnt_s}{cnt_f});
        end
        for cnt_iners = 1:length(VAL_iner{:,1})
            temp = [];
            for cnt_s = 1:N
                if size(VAL_iner{cnt_s}{cnt_iners},1) < size(VAL_iner{cnt_s}{cnt_iners},2)
                    temp = [temp, VAL_iner{cnt_s}{cnt_iners}];
                else
                    temp = [temp, VAL_iner{cnt_s}{cnt_iners}];
                end
            end
            Y.(FN{1}{cnt_f}).(FN_iner{:,cnt_s}{cnt_iners}) = temp;
        end
        clear VAL_iner FN_iner;
    end
end
end

function [c] = appendStruct(a,b)
%APPENDSTRUCT Appends two structures ignoring duplicates
%   Developed to append two structs while handling cases of non-unique
%   fieldnames.  The default keeps the last occurance of the duplicates in
%   the appended structure.
    ab = [struct2cell(a); struct2cell(b)];
    abNames = [fieldnames(a); fieldnames(b)];
    [~,iab] = unique(abNames,'last');
    c = cell2struct(ab(iab),abNames(iab));
end
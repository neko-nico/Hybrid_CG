%保存matlab脚本文件副本的代码
%定义输出文件夹
folderName = 'Test';

% 如果文件夹不存在，则创建
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% 获取当前时间（格式：YYYYmmdd-HHMMSS）
time_format = datetime('now','Format','MMdd-HHmmss');
timeStr = char(time_format);   % 转成字符串

% 拼接日志文件名
logFileName = fullfile(folderName, [timeStr '_log' '.txt']);
% 打开文件（写入模式，若文件已存在则覆盖）
logfid = fopen(logFileName,'w');

% 判断文件是否成功打开
if logfid == -1
    error('无法创建日志文件');
end

% 保存程序文件副本
% 获取当前脚本完整路径
scriptName = [mfilename(), '.m'];
% 复制并重命名
copyfile(scriptName, fullfile(folderName, [timeStr, '_', scriptName]));

fprintf(logfid, '副本已保存\n');
fclose(logfid);
clear,clc;
close all;  % 安全关闭所有绘图窗口
%opt-----------------------------------------
WGS84_a = 6378137.0;%椭球长半轴
WGS84_f = 1.0/298.257223563;%椭球扁率
WGS84_e2 = 0.00669437999014;%椭球第一偏心率的平方
Trans_Switch = 0;% 0:BLH->XYZ / 1:XYZ->BLH
%trans-----------------------------------------
Pre_B='30 31 54.0';
Pre_L='114 21 26.2';
Pre_H=28.2;
Pre_X=0;
Pre_Y=0;
Pre_Z=0;
B=0.0;
L=0.0;
H=0.0;
X=0.0;
Y=0.0;
Z=0.0;
if(Trans_Switch)
    [B,L,H]=blh2xyz(Pre_X,Pre_Y,Pre_Z);
    write_result(B,L,H);
else
    Rad_B=dms2rad(Pre_B);
    Rad_L=dms2rad(Pre_L);
    [X,Y,Z]=blh2xyz(Rad_B,Rad_L,Pre_H,WGS84_a,WGS84_e2);
    write_result(X,Y,Z);
end
function [X, Y, Z] = blh2xyz(B, L, H, a, e2)
% 输入：B, L（弧度），H（高度），a（长半轴），e2（第一偏心率平方）
    N = a / sqrt(1 - e2 * sin(B) * sin(B));
    X = (N + H) * cos(B) * cos(L);
    Y = (N + H) * cos(B) * sin(L);
    Z = (N*(1 - e2) + H) * sin(B);
end

function [B, L, H] = xyz2blh(X, Y, Z, a, e2)
% 使用 Bowring 法快速反解 BLH
    p = sqrt(X.^2 + Y.^2);
    r = sqrt(p.^2 + Z.^2);
    tan_u = (Z / p) * (a / sqrt(a^2 - e2 * Z.^2));
    tan_B = (Z + e2 * a / sqrt(1 - e2 * sin(tan_u).^2) * sin(tan_u)) / p;
    B = atan(tan_B);
    L = atan2(Y, X);
    N = a ./ sqrt(1 - e2 * sin(B).^2);
    H = p ./ cos(B) - N;
end

function rad = dms2rad(angleStr)
% DMS2RAD 将 '度 分 秒' 格式的字符串转换为弧度值
% 输入:
%   angleStr - 字符串，格式如 '30 31 54.0'
%
% 输出:
%   rad - 对应的角度值（以弧度表示）

    % 拆分字符串
    parts = sscanf(angleStr, '%f%f%f');
    
    % 提取度、分、秒
    degrees = parts(1);
    minutes = parts(2);
    seconds = parts(3);
    
    % 转换为十进制度
    decimalDegrees = degrees*1.0 + minutes/60.0 + seconds/3600.0;
    
    % 转换为弧度
    rad = deg2rad(decimalDegrees);
end

function write_result(X, Y, Z)
% write_result - 将三个数字写入 result.txt，每个数字一行，保留四位小数，不使用科学计数法
%
% 用法：
%   write_result(X, Y, Z)
%
% 示例：
%   write_result(123456.78901, 0.1234, -56.78)

    % 打开文件用于写入（覆盖旧内容）
    fid = fopen('Trans-result.txt', 'w');

    if fid == -1
        error('无法打开 result.txt');
    end

    % 写入每个数字，保留4位小数，不使用科学计数法
    fprintf(fid, '%.4f\n', X);
    fprintf(fid, '%.4f\n', Y);
    fprintf(fid, '%.4f\n', Z);

    fclose(fid);
end

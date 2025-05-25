@echo off
setlocal

:: 设置要删除的文件类型
set "ext1=.trace"
set "ext2=.pospsu_res"
set "ext3=.pospsu_snr"
set "ext4=.pos"

:: 获取当前路径
set "target_dir=%cd%"

:: 切换到目标路径
cd /d "%target_dir%"

echo 正在删除以下类型的文件：
echo   %ext1%
echo   %ext2%
echo   %ext3%
echo.

:: 删除 .trace 文件
if exist *%ext1% (
    echo 删除 %ext1% 文件...
    del /q *%ext1%
) else (
    echo 没有找到 %ext1% 文件。
)

:: 删除 .pospsu_res 文件
if exist *%ext2% (
    echo 删除 %ext2% 文件...
    del /q *%ext2%
) else (
    echo 没有找到 %ext2% 文件。
)

:: 删除 .pospsu_snr 文件
if exist *%ext3% (
    echo 删除 %ext3% 文件...
    del /q *%ext3%
) else (
    echo 没有找到 %ext3% 文件。
)

echo.
echo 所有指定格式文件已清理完成！
echo.
:: 将结果pos文件转移至Solution文件夹
:: 创建 Solution 文件夹（如果不存在）
if not exist "Solution" (
    mkdir "Solution"
)
:: 将所有结果 .pos 文件移动到 Solution 文件夹
echo 正在转移 %ext4% 文件。
if exist *%ext4% (
    move *.pos "Solution\" >nul 2>&1
) else (
    echo 没有找到 %ext4% 文件。
)

echo.
echo 全部.pos文件均已移动到 Solution 文件夹
:: 保持窗口不关闭
cmd /k
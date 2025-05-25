@echo off
setlocal

:: 获取当前目录
set "target_dir=%cd%"

:: 切换到目标目录
cd /d "%target_dir%"

:: 创建 RawFile 文件夹（如果不存在）
if not exist "RawFile" (
    mkdir "RawFile"
)

:: 将 .crx 文件重命名为 .25d
ren *.crx *.25d

:: 遍历所有 .25d 文件并执行 crx2rnx 命令
for %%f in (*.25d) do (
    echo 正在处理文件：%%f
    crx2rnx "%%f"
)

:: 将所有原始 .25d(原.crx) 文件移动到 RawFile 文件夹
move *.25d "RawFile\" >nul 2>&1

echo.
echo 所有文件处理完成，原始文件已移动至 RawFile 文件夹。
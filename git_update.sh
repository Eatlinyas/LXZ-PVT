#!/bin/bash

# ================ 配置区域 =================
# 要提交的文件夹路径（相对于项目根目录）
FOLDER="src"

# 提交信息（可根据需要修改）
COMMIT_MSG="更新 $FOLDER 内容"

# Git 远程仓库分支名（通常是 main 或 master）
BRANCH="main"
# ===========================================

# 检查是否在 Git 项目中
if [ ! -d ".git" ]; then
  echo "❌ 错误：当前目录不是一个 Git 仓库。"
  exit 1
fi

# 检查目标文件夹是否存在
if [ ! -d "$FOLDER" ]; then
  echo "❌ 错误：目标文件夹 '$FOLDER' 不存在。"
  exit 1
fi

echo "✅ 正在提交 '$FOLDER' 文件夹的修改..."

# 添加指定文件夹的改动
git add "$FOLDER/"

# 检查是否有改动被添加
CHANGES=$(git diff --cached --name-only)
if [ -z "$CHANGES" ]; then
  echo "⚠️ 没有检测到 '$FOLDER' 中的改动，无需提交。"
  exit 0
fi

# 提交改动
git commit -m "$COMMIT_MSG"

# 推送到远程仓库
git push origin "$BRANCH"

echo "🎉 成功将 '$FOLDER' 的修改提交并推送到远程仓库 '$BRANCH' 分支。"

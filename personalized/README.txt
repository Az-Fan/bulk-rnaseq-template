个性化分析模块
==============

用途
----
保存特定绘图、补充统计、临时假设检验或项目专属分析。

目录
----
code/       每个分析一个独立 R 脚本
work/       中间文件
results/    最终结果

规则
----
1. 不修改 scripts/pipeline 中的标准 Step 01–06。
2. 不把个性化结果复制到标准 results/<ProjectName>。
3. 不进入标准 manifest.csv。
4. 每个分析使用唯一名称，例如 endothelial_signature_plot。
5. 每个结果目录必须有 README.txt，记录输入、方法和阈值。

运行示例
--------
.\run_personalized.ps1 -Name example_analysis -Config config/project.psd1

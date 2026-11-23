# 中文使用说明

如果要核对论文图表，建议从 S5 冻结分析开始。若要检查某一步对结果的影响，再从原始 reads、MAG 或模型输入重新运行对应的分支。公开仓库提供代码和使用说明，测序数据从 NCBI 取得，绘图表使用论文附件或作者本地保存的输入。

**从附件重绘。** 将 S1–S4 四个工作簿、`Supplementary_Data_S3_Fig3_Addendum.zip`、`Supplementary_Data_S5.zip` 和 `Supplementary_Data_S6_FCS.zip` 放入仓库同级的 `books/`。在仓库目录中执行下列命令。Windows PowerShell 可直接逐行运行。Rscript 需要在 PATH 中。

```text
conda env create -f env/replay.yml
conda activate zabuye-replay
python -m unittest discover -s tests -v
python workflow/stage.py --mode replay --s5 ../books/Supplementary_Data_S5.zip --s6 ../books/Supplementary_Data_S6_FCS.zip --project ../replay
python workflow/paper_inputs.py --books ../books --project ../replay --out ../paper-inputs
python workflow/figures.py --group paper --books ../books --project ../replay --local-inputs ../paper-inputs --out ../paper-figures
```

这条路线可重绘 5 张主图中的 4 张，以及 28 张补图中的 27 张。Figure 4 和 Figure S9 各需一张保留在作者本地的蛋白 pI 密度表，两者分别对应全体与高质量 Zabuye 水体 MAG。具备这两张表时运行 `--group all` 即可得到全部 33 个图号。只核对 S20、S21 等修回新增图时，可单独运行对应组。Figure 3C 保留全部 12 列。[逐图输入表](figures.md)

**从原始测序重跑。** 从 BioProject 导出 RunInfo，使用其中的 `LibraryName` 与 `Run`，再由 `manifest.py` 对照 S5 的 MEGAHIT 输入顺序生成 35 个文库的清单。计算流程包含 reads 质控、分类和群落 KO、组装分箱、MAG 质控和注释、覆盖度、配对模型、系统树及绘图。首次运行 `run.py` 默认只列命令，加入 `--execute` 才执行。原始四组共组装使用独立保存的预清洗 FASTQ，当前缺少证据证明其与公开 reads 经 fastp 得到的文件完全相同。对应的历史分析、同输入 14 文库敏感性和冻结结果复算，在运行指南中给出了各自的输入。[NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1531265) · [Linux 与超算运行指南](compute.md)

**比较结果时保留分析单位。** 日期比较限于 2023 年 4 月与 10 月的七个配对站点。SCA、GMP 和检出携带者比例使用不同分母。MAG 功能分解描述加性归属，蛋白 pI 是序列性状。外部沉积物数据的环境模型以湖为单位。手工制图部分可通过 `artwork.py` 从自有 Word 稿件导出终稿图像，与代码重绘结果对照。

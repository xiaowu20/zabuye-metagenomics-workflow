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

这条路线可重绘 5 张主图中的 4 张，以及 28 张补图中的 27 张。Figure 4 和 Figure S9 还需要全体与高质量 Zabuye 水体 MAG 的等权蛋白 pI 密度曲线。可直接使用本地归档的曲线，也可用公开的 `workflow/protein_density.R` 从逐蛋白 pI 和 S4 重建，再执行 `--group all` 生成全部 33 个图号。只核对 S20、S21 等修回新增图时，可单独运行对应组。Figure 3C 保留全部 12 列。[逐图输入与密度重算命令](figures.md)

**从原始测序重跑。** 从 BioProject 导出 RunInfo，用 `manifest.py` 将 `LibraryName` 和 `Run` 与 S5 中的 MEGAHIT 日志和参数文件对应起来。新清单记录 35 个文库各自的输入顺序和双端合计读数。`run.py` 默认打印命令，加入 `--execute` 才执行。`reads` 分支用 fastp 1.3.2 处理用于分类、功能分析、分箱深度估计及新增 14 文库组装的 reads，未启用重叠区纠错。[NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1531265)

四组历史共组装另用 fastp 0.22.0，开启双端接头检测和重叠区纠错，合格碱基阈值为 Phred 20，最短读长为 50 bp。4 月 9 日的 D1 日志记录了随后 MEGAHIT 读取的两个输出路径。按上述参数用 8 线程重新处理作者归档的 35 个原始文库，各文库过滤后的读数均与历史组装日志一致。历史 FASTQ 缺少可供比对的内容校验值。这次核验覆盖原始文件到质控读数，未从 SRA 重新导出文件并验证，也未重跑后续全部组装、分箱和绘图。

执行 `reads` 后，`historical --historical-fastp /path/to/fastp-0.22.0` 默认从下载的原始双端文件重建历史输入。要用作者归档的原始 FASTQ，可加 `--historical-raw DIR`。新清单会在组装前逐库核对过滤读数。已有 `00_qc_raw` 过滤后文件时，可直接使用 `historical --assembly-reads DIR`。[完整命令、文件名及环境设置](compute.md#executing-the-raw-data-branches)

**比较结果时保留分析单位。** 日期比较限于 2023 年 4 月与 10 月的七个配对站点。SCA、GMP 和检出携带者比例使用不同分母。MAG 功能分解描述加性归属，蛋白 pI 是序列性状。外部沉积物数据的环境模型以湖为单位。手工制图部分可通过 `artwork.py` 从自有 Word 稿件导出终稿图像，与代码重绘结果对照。

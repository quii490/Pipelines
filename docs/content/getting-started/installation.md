# 环境与检查

统一激活MCC目录下的conda环境。
```bash
## RNA-seq pipeline
conda activate /home/machicheng/.conda/envs/rnaseq
```
```bash
## ATAC-seq和CUT&RUN pipeline
conda activate /home/machicheng/.conda/envs/chipseq
```
### 最小入口检查（可选）

```bash
bash RNA-seq/rnaseq/run_auto_rnaseq.sh --help
bash ATAC-seq/run_auto_atacseq.sh --help
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
# 这只证明入口可读取，不证明所有依赖已安装。正式运行前还应检查：
java -version
nextflow -version
conda --version
samtools --version
```
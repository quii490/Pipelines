# 命令行最小知识
掌握下面这些操作即可安全使用流程。
## 路径与文件

```bash
pwd                         # 当前目录
cd /path/to/project #进入目录
ls -lh /path/to/project     # 查看文件和大小
find /path/to/fastq -type f -name '*.fastq.gz'
test -r /path/to/file && echo readable
test -w /path/to/results && echo writable
df -h      # 可用磁盘

mv /path/to/project1 /path/to/project2 #移动文件
cp /path/to/project1 /path/to/project2 #复制文件
rm /path/to/project #删除文件。慎用rm -rf!!!!
```

有空格的路径必须整体加引号。正式分析更建议避免空格和中文标点。
## 运行代码
```
bash脚本，bash /path/to/project.sh
python脚本，python /path/to/project.py
```
## 变量让命令更清楚
也可以不用。平时用不到变量。
```bash
PROJECT=/path/to/project
FASTQ="$PROJECT/fastq"
RESULTS="$PROJECT/results"
`echo "$RESULTS"` #检查变量
```
### glob 要加引号

```bash
--bw-glob "/path/to/bw/*.bw"
```

引号让脚本自己处理 glob。若不加引号，shell 可能提前展开成多个参数，导致“unknown option”或只读取第一个文件。



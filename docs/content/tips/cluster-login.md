# 登录和使用计算集群

## 第一次登录

```bash
ssh user@login.example.org
```

首次连接会显示 host key fingerprint。不要看到提示就直接输入 `yes`；先通过实验室管理员或可信内部渠道核对 fingerprint。密码不会在终端显示字符，这是正常现象。

## 推荐使用 SSH key

在自己的电脑生成 key，不要在共享服务器生成后把 private key 下载回来：

```bash
ssh-keygen -t ed25519 -a 100 -C "your-name@lab"
ssh-copy-id user@login.example.org
```

macOS 没有 `ssh-copy-id` 时，把 `~/.ssh/id_ed25519.pub` 的**公钥**内容交给管理员或按机构流程添加。绝不发送 `id_ed25519` 私钥。

为 key 设置 passphrase。若使用 agent：

```bash
ssh-add ~/.ssh/id_ed25519
ssh-add -l
```

## 用 SSH config 简化命令

本机 `~/.ssh/config`：

```sshconfig
Host lab-cluster
    HostName login.example.org
    User user
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 60
    ServerAliveCountMax 3
```

以后：

```bash
ssh lab-cluster
```

该文件可能含内部主机名，不要提交到公开 GitHub。共享模板使用 `config.example` 和假地址。

## 登录节点不能做什么

Login node 用于：编辑文本、提交作业、查看状态、传小文件。不要在 login node 上直接运行 STAR、Bowtie2、samtools sort、deepTools、R 大矩阵或完整 pipeline。

需要交互调试时申请 compute node，例如 Slurm：

```bash
srun --partition=compute --cpus-per-task=4 --mem=16G --time=01:00:00 --pty bash
hostname
```

分区和限制由集群决定，以上名称只是示例。

## 跳板机

若机构要求 bastion，可在 SSH config 使用：

```sshconfig
Host lab-cluster
    HostName compute-login.internal
    User user
    ProxyJump user@bastion.example.org
```

不要用文档、AI 对话或 GitHub issue 传播真实内网拓扑。

## 连接失败快速检查

```bash
ssh -v lab-cluster
```

只分享经过脱敏的末尾日志。常见原因：VPN 未连接、账户过期、key 权限过宽、host key 变化、网络限制或用户名错误。host key 突然变化时先联系管理员，不能直接删除 `known_hosts` 记录来绕过警告。

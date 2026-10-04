# Lab 7：测试、调试与性能分析

课程：C 程序与算法分析项目实践课  
作者：Gong Helin  
年份：2026

本代码包用于第七节课。程序围绕矩阵求解库的正确性验证、错误路径检查和性能测试组织，覆盖矩阵基本运算、矩阵乘法、高斯消元、LU 分解和线性方程组求解。

## 1. 文件结构

```text
lec7code/
  matrix_core.h / matrix_core.c       矩阵数据结构、内存管理和通用函数
  matrix_ops.h  / matrix_ops.c        矩阵基本运算、矩阵乘法和范数
  matrix_solve.h / matrix_solve.c     Gaussian elimination 线性方程组求解
  matrix_lu.h   / matrix_lu.c         LU 分解、前代、回代和行列式
  test_runner.c                      单元测试和错误路径测试
  benchmark.c                        性能测试：矩阵乘法、高斯消元和 LU 分解
  student_practice.c                 课内练习程序
  Makefile                           编译和运行规则
```

## 2. 编译方法

在终端进入代码目录，执行：

```bash
make
```

该命令生成三个可执行程序：

```text
lab7_tests       单元测试程序，由 test_runner.c 编译得到
lab7_benchmark   性能测试程序，由 benchmark.c 编译得到
lab7_practice    课内练习程序，由 student_practice.c 编译得到
```

## 3. 运行方法

运行单元测试和错误路径测试：

```bash
make run-tests
```

运行性能测试：

```bash
make run-bench
```

运行课内练习程序：

```bash
make run-practice
```

清理编译生成文件：

```bash
make clean
```

## 4. 测试内容

`test_runner.c` 按以下内容组织测试：

1. 矩阵加法手算样例验证；
2. 矩阵乘法手算样例验证；
3. 维度不匹配和下标越界错误码检查；
4. 行交换函数 `MatrixSwapRows` 检查；
5. 无主元高斯消元求解检查；
6. 列主元高斯消元求解检查；
7. 奇异矩阵在高斯消元中的错误码检查；
8. LU 分解后用 `L * U` 重构原矩阵；
9. 利用 LU 分解求解 `Ax = b`；
10. 奇异矩阵在 LU 分解中的错误码检查。

## 5. 性能测试内容

`benchmark.c` 是性能测试程序，执行命令为：

```bash
make run-bench
```

该程序包含两个测试入口：

```c
BenchmarkMatrixMultiplication();
BenchmarkLinearSolvers();
```

其中：

- `BenchmarkMatrixMultiplication` 比较 `MatrixMultiply` 和 `MatrixMultiplyIKJ` 的平均运行时间；
- `BenchmarkLinearSolvers` 比较 `GaussianSolvePartialPivot` 和 `LUDecomposeNoPivot` 的平均运行时间；
- 每个规模重复运行 3 次，输出平均时间。

典型输出格式如下：

```text
[Benchmark 1] Matrix multiplication
repeats = 3
n =  100 | ijk = ... s | ikj = ... s

[Benchmark 2] Linear solver kernels
repeats = 3
n =   50 | Gaussian(partial pivot) = ... s | LU = ... s
```

性能数据与电脑配置、编译选项、系统负载和矩阵规模有关。记录报告数据时，应同时写明测试环境、矩阵规模、重复次数和平均时间。

## 6. 推荐使用顺序

```bash
make clean
make
make run-tests
make run-bench
make run-practice
make clean
```

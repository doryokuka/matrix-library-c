# =============================================================================
# File   : Makefile
# Course : C Program and Algorithm Analysis Lab -- Matrix Solver Library
# Author : Gong Helin
# Year   : 2026
# Note   : 第七节课构建文件：单元测试、错误路径测试和性能测试。
# =============================================================================

CC = gcc

# C99 标准；打开常用警告；O0 和 -g 便于调试与定位问题。
CFLAGS = -std=c99 -Wall -Wextra -pedantic -O0 -g

# fabs 等数学函数位于 libm。
LDFLAGS = -lm

TEST_TARGET = lab7_tests
BENCH_TARGET = lab7_benchmark
PRACTICE_TARGET = lab7_practice

# 矩阵库公共模块：核心数据结构、基本运算、高斯消元、LU 分解和 Jacobi 对角化。
COMMON_OBJS = matrix_core.o matrix_ops.o matrix_solve.o matrix_lu.o matrix_diag.o

all: $(TEST_TARGET) $(BENCH_TARGET) $(PRACTICE_TARGET)

$(TEST_TARGET): test_runner.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_TARGET) test_runner.o $(COMMON_OBJS) $(LDFLAGS)

$(BENCH_TARGET): benchmark.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) -o $(BENCH_TARGET) benchmark.o $(COMMON_OBJS) $(LDFLAGS)

$(PRACTICE_TARGET): student_practice.o $(COMMON_OBJS)
	$(CC) $(CFLAGS) -o $(PRACTICE_TARGET) student_practice.o $(COMMON_OBJS) $(LDFLAGS)

test_runner.o: test_runner.c matrix_core.h matrix_ops.h matrix_solve.h matrix_lu.h matrix_diag.h
	$(CC) $(CFLAGS) -c test_runner.c

benchmark.o: benchmark.c matrix_core.h matrix_ops.h matrix_solve.h matrix_lu.h matrix_diag.h
	$(CC) $(CFLAGS) -c benchmark.c

student_practice.o: student_practice.c matrix_core.h matrix_ops.h matrix_solve.h matrix_lu.h matrix_diag.h
	$(CC) $(CFLAGS) -c student_practice.c

matrix_core.o: matrix_core.c matrix_core.h
	$(CC) $(CFLAGS) -c matrix_core.c

matrix_ops.o: matrix_ops.c matrix_ops.h matrix_core.h
	$(CC) $(CFLAGS) -c matrix_ops.c

matrix_solve.o: matrix_solve.c matrix_solve.h matrix_core.h
	$(CC) $(CFLAGS) -c matrix_solve.c

matrix_lu.o: matrix_lu.c matrix_lu.h matrix_core.h
	$(CC) $(CFLAGS) -c matrix_lu.c

matrix_diag.o: matrix_diag.c matrix_diag.h matrix_core.h matrix_ops.h
	$(CC) $(CFLAGS) -c matrix_diag.c

run-tests: $(TEST_TARGET)
	./$(TEST_TARGET)

run-bench: $(BENCH_TARGET)
	./$(BENCH_TARGET)

run-practice: $(PRACTICE_TARGET)
	./$(PRACTICE_TARGET)

clean:
	rm -f $(TEST_TARGET) $(BENCH_TARGET) $(PRACTICE_TARGET) *.o

#!/bin/bash

BUILD_DIR="build_release"
OUTPUT_NAME_PRE="benchmark_wo_semijoin_reduction_k_"
OUTPUT_NAME_POST=".json"

KS=(10 50 100 1000 5000)

for K in KS
do
	echo "Running Benchmark For K = ${K}"
	./${BUILD_DIR}/hyriseBenchmarkJoinOrder -r 100 -o "${OUTPUT_NAME_PRE}${K}${OUTPUT_NAME_POST}" --top_k ${K}
done


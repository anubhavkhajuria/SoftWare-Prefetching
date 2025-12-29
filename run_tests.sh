rm testcases/pagerank_vector.ll
rm testcases/pagerank_vector_preprocess.ll
rm testcases/pagerank_vector_prefetch_out.ll
rm testcases/pagerank_vector_opt.o
rm testcases/pagerank_vector_exe

rm testcases/pagerank_array.ll
rm testcases/pagerank_array_preprocess.ll
rm testcases/pagerank_array_prefetch_out.ll
rm testcases/pagerank_array_opt.o
rm testcases/pagerank_array_exe

clang++ -emit-llvm -S testcases/pagerank_vector.cpp -o testcases/pagerank_vector.ll -std=c++17 -Xclang -disable-O0-optnone -O0 -fno-inline
opt -passes='mem2reg,loop-simplify,lcssa,indvars,simplifycfg' -S testcases/pagerank_vector.ll -o testcases/pagerank_vector_preprocess.ll
opt -load-pass-plugin ./build/pass/prefetchIMAPass.dylib -passes='prefetch-ima' -ima-prefetch-distance=16 -ima-degree-threshold=16 -S testcases/pagerank_vector_preprocess.ll -o testcases/pagerank_vector_prefetch_out.ll
clang++ -O3 -c testcases/pagerank_vector_prefetch_out.ll -o testcases/pagerank_vector_opt.o
clang++ testcases/pagerank_vector_opt.o -O3 -o testcases/pagerank_vector_exe
./testcases/pagerank_vector_exe Datasets/web-Stanford.txt 0 9000000

clang++ -emit-llvm -S testcases/pagerank_array.cpp -o testcases/pagerank_array.ll -std=c++17 -Xclang -disable-O0-optnone -O0 -fno-inline
opt -passes='mem2reg,loop-simplify,lcssa,indvars,simplifycfg' -S testcases/pagerank_array.ll -o testcases/pagerank_array_preprocess.ll
opt -load-pass-plugin ./build/pass/prefetchIMAPass.dylib -passes='prefetch-ima' -ima-prefetch-distance=16 -ima-degree-threshold=16 -S testcases/pagerank_array_preprocess.ll -o testcases/pagerank_array_prefetch_out.ll
clang++ -O3 -c testcases/pagerank_array_prefetch_out.ll -o testcases/pagerank_array_opt.o
clang++ testcases/pagerank_array_opt.o -O3 -o testcases/pagerank_array_exe
./testcases/pagerank_array_exe Datasets/web-Stanford.txt 0 9000000

# For linux based machine, replace prefetchIMAPass.so with prefetchIMAPass.dylib
# For mac based machine, replace prefetchIMAPass.dylib with prefetchIMAPass.so

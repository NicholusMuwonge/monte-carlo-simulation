If you're running this on Mac M1 chip and above;

`arch -x86_64 /usr/local/opt/llvm/bin/clang -Xclang -fopenmp parallel_programming.c -o parallel_programming -L/usr/local/opt/llvm/lib -lomp`# monte-carlo-simulation

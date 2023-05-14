if [ $# -eq 0 ]
  then
    echo "[!] Please specify input file number (e.g.: ./run.sh 0)"
else
    echo "[*] Compiling...";
    g++ -O3 ../src/main.cpp ../src/sssp.cpp ../src/utils.cpp ../src/graph.cpp -o SSSP && 
    echo "[*] Running..." && 
    (time ./SSSP < input/input$1.txt 2>debug.out && echo "[*] Done" || echo "Input file named 'input$1.txt' not found in input folder.")
fi

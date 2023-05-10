if [ $# -eq 0 ]
  then
    echo "Please specify input file number (e.g.: ./run.sh 0)"
else
    g++ ../src/main.cpp ../src/sssp.cpp ../src/utils.cpp ../src/graph.cpp -o SSSP && ./SSSP < input$1.txt 2>debug.out || echo "Input file named 'input$1.txt' not found here."
fi

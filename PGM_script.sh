#!/bin/bash

# Function to run the program.
run_program() {
    local d="$1"
    local L="$2"

    echo "Running program with d=${d} and L=${L}..."

    gcc -std=c11 -ferror-limit=0 main_PGM.c func.c mt19937ar.c -I/opt/homebrew/opt/openblas/include/ -L/opt/homebrew/opt/openblas/lib -Wl,-rpath,/opt/homebrew/opt/openblas/lib -lopenblas -lm

    ./a.out "$d" "$L"
}

export -f run_program

# Values for d.
d_values=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22 23 24 26 27 28 29 30)

# Values for L.
L_values=(1000)

# Check for parallel flag.
if [[ $1 == "--parallel" ]]; then
    # Run the program in parallel.
    parallel run_program ::: "${d_values[@]}" ::: "${L_values[@]}"
else
    # Run the program sequentially.
    for d in "${d_values[@]}"; do
        for L in "${L_values[@]}"; do
            run_program "$d" "$L"
        done
    done
fi

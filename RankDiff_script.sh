#!/bin/bash

# Define the function for running the program
run_program() {
  local d="$1"
  local L="$2"
  local N="$3"
  local I="$4"
  local r="$5"  # Add the new variable 'r'

  # Generate a unique name for each executable
  local unique_identifier="d${d}_L${L}_N${N}_I${I}_r${r}"  # Add 'r' here
  local output_executable="P_${unique_identifier}"

  echo "Running program with d=${d}, L=${L}, N=${N}, I=${I}, r=${r}, Executable=${output_executable}"

  # Compile the code and specify a unique output executable
  # For macOS (commented out, enable if you are running on macOS)
  gcc -std=c11 main_RMB_RankDiff.c func.c mt19937ar.c -I/opt/homebrew/opt/openblas/include/ -L/opt/homebrew/opt/openblas/lib -Wl,-rpath,/opt/homebrew/opt/openblas/lib -lopenblas -lm -o "${output_executable}"

  # For Linux
  #gcc -std=c11 main_RMB_RankDiff.c func.c mt19937ar.c -Wl,-z,muldefs -llapacke -llapack -lopenblas -lm -o "${output_executable}"

  # Run the program
  "./${output_executable}" "$d" "$L" "$N" "$I" "$r"  # Add 'r' here

  # Delete the unique executable after the job is done
  rm -f "${output_executable}"
}
export -f run_program

# Values for d
d_values=(2 4 8 16)

# Values for L
L_values=(100000)

# Values for N
N_values=(1 10 100)

# Values for I
I_values=($(seq 1 1 1000))

# Check for parallel flag
if [[ $1 == "--parallel" ]]; then
  # Run the program in parallel
  for d in "${d_values[@]}"; do
    r_values=$(seq 1 $d)
    parallel --jobs "$B_JOBS" run_program ::: "$d" ::: "${L_values[@]}" ::: "${N_values[@]}" ::: "${I_values[@]}" ::: "${r_values[@]}"
  done
else
  # Run the program sequentially
  for I in "${I_values[@]}"; do
    for d in "${d_values[@]}"; do
      for L in "${L_values[@]}"; do
        for N in "${N_values[@]}"; do
          for r in $(seq 1 $d); do  # Looping through 'r' values in the range [1,d]
            run_program "$d" "$L" "$N" "$I" "$r"  # Add 'r' here
          done
        done
      done
    done
  done
fi

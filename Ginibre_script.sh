#!/bin/bash

# Function to run the program.
run_program() {
  local d="$1"
  local L="$2"
  local N="$3"
  local I="$4"

  # Generate a unique name for each executable
  local unique_identifier="d${d}_L${L}_N${N}_I${I}"
  local output_executable="P_${unique_identifier}"

  echo "Running program with d=${d}, L=${L}, N=${N} and I=${I}, Executable=${output_executable}"

  # Compile the code and specify a unique output executable
  # For macOS
  #gcc -std=c11 main_RMB_Ginibre.c func.c mt19937ar.c -I/opt/homebrew/opt/openblas/include/ -L/opt/homebrew/opt/openblas/lib -Wl,-rpath,/opt/homebrew/opt/openblas/lib -lopenblas -lm -o "${output_executable}"

  # For Linux
  gcc -std=c11 main_RMB_Ginibre.c func.c mt19937ar.c -Wl,-z,muldefs -llapacke -llapack -lopenblas -lm -o "${output_executable}"

  # Run the program
  "./${output_executable}" "$d" "$L" "$N" "$I"

  # Delete the unique executable after the job is done
  rm -f "${output_executable}"
}
export -f run_program

# Values for d.
d_values=(2 4 8 16)
#d_values=(2)

# Values for L.
L_values=(100000)

# Value for N.
N_values=(1 10 100)
#N_values=(100)

# Value for I
I_values=($(seq 1 1 1000))
#I_values=(58 87 95 34 55 61 66 73 85)
#I_values=(1 7 17 24 41 57 79 84 93 97 10 18 19 36 38 49 55 59 103 118 124 290 334 394 396 403)

# Check for parallel flag.
if [[ $1 == "--parallel" ]]; then
  # Run the program in parallel.
  parallel --jobs "$B_JOBS" run_program ::: "${d_values[@]}" ::: "${L_values[@]}" ::: "${N_values[@]}" ::: "${I_values[@]}"
else
  # Run the program sequentially.
  for I in "${I_values[@]}"; do
    for d in "${d_values[@]}"; do
      for L in "${L_values[@]}"; do
        for N in "${N_values[@]}"; do
          run_program "$d" "$L" "$N" "$I"
        done
      done
    done
  done
fi

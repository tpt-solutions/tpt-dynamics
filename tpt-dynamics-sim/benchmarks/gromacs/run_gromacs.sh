#!/bin/bash
cd ~/gmx-bench || exit 1
export LD_LIBRARY_PATH=~/scratch-gmx/gmx/lib:$LD_LIBRARY_PATH
GMX=~/scratch-gmx/gmx/bin/gmx
for n in 1 4 8 16 20; do
  $GMX mdrun -s lj.tpr -ntmpi 1 -ntomp "$n" -nsteps 2000 -noconfout -nb cpu -pin on -g "pin-$n.log" -e "pin-$n.edr" > "pin-$n.out" 2>&1
  echo "nt=$n exit=$?"
  grep -E "Performance" -A2 "pin-$n.out" | head -4
done

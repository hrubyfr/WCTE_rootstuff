
COUNTER=".run_count"

if [ ! -f "$COUNTER" ]; then
	echo 0 > "$COUNTER"
fi

count=$(cat "$COUNTER")

count=$((count+1))

echo $count > "$COUNTER"

OUTPUT="increasing_energy_${count}.dat"

declare -A energies
runs=()

input_file="${1:?Usage: $0 <run_energies_file>}"

while read -r run energy; do
    energies[$run]=$energy
    runs+=($run)
done < "$input_file"

for run in "${runs[@]}"; do
	echo "Processing file $run"


	energy="${energies[$run]}"
	file=$(ls /media/frantisek/T7\ Shield/WCTE_data/WCTE_offline_R${run}S0_*.root 2>/dev/null)

	echo "Energy of the run $run is $energy"
	echo "running filename $file"

	if [ -z "$file" ]; then
		echo "⚠️  Warning: File for run $run not found"
		continue
	fi
	root -b -q "read_matched_data.cc+(\"$file\", \"$OUTPUT\", $energy)"
done

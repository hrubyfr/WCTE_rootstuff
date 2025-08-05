
COUNTER=".run_count"

if [ ! -f "$COUNTER" ]; then
	echo 0 > "$COUNTER"
fi

count=$(cat "$COUNTER")

count=$((count+1))

echo $count > "$COUNTER"

OUTPUT="measured_sigma_${count}.dat"

declare -A energies
runs=()

while read -r run energy; do
    energies[$run]=$energy
    runs+=($run)
done < run_energies.dat

files=("WCTE_offline_R1362S0_VME1443.root" "WCTE_offline_R1373S0_VME1447.root" "WCTE_offline_R1374S0_VME1448.root" "WCTE_offline_R1375S0_VME1450.root" "WCTE_offline_R1376S0_VME1451.root" "WCTE_offline_R1398S0_VME1456.root" "WCTE_offline_R1518S0_VME1532.root" "WCTE_offline_R1521S0_VME1534.root" "WCTE_offline_R1541S0_VME1544.root" "WCTE_offline_R1556S0_VME1547.root" "WCTE_offline_R1536S0_VME1540.root" "WCTE_offline_R1537S0_VME1541.root" "WCTE_offline_R2158S0_VME1879.root" "WCTE_offline_R1510S0_VME1529.root" "WCTE_offline_R1514S0_VME1531.root")

for file in "${files[@]}"; do
	echo "Processing file $file"

	run_part="${file#*R}"
	run_number="${run_part%%S*}" 

	energy="${energies[$run_number]}"

	echo $energy

	if [ -z "$energy" ]; then
		echo "Warning: No energy found for run $run_number! Skipping file $file"
		continue
	fi

	root -b -q "read_matched_data.cc+(\"/media/frantisek/T7 Shield/WCTE_data/$file\", \"$OUTPUT\", \"$energy\")"
done

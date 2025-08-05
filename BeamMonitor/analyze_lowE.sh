
COUNTER=".run_count"

if [ ! -f "$COUNTER" ]; then
	echo 0 > "$COUNTER"
fi

count=$(cat "$COUNTER")

count=$((count+1))

echo $count > "$COUNTER"

OUTPUT="lowE_sigma_${count}.dat"

files=( "WCTE_offline_R1541S0_VME1544.root" "WCTE_offline_R1536S0_VME1540.root" "WCTE_offline_R1537S0_VME1541.root" "WCTE_offline_R2158S0_VME1879.root" "WCTE_offline_R1510S0_VME1529.root" "WCTE_offline_R1582S0_VME1562.root" "WCTE_offline_R1586S0_VME1563.root" "WCTE_offline_R1893S0_VME1741.root" "WCTE_offline_R2161S0_VME1881.root" "WCTE_offline_R2163S0_VME1882.root" "WCTE_offline_R2165S0_VME1883.root")

for file in "${files[@]}"; do
	echo "Processing file $file"
	root -b -q "read_matched_data.cc+(\"/media/frantisek/T7 Shield/WCTE_data/$file\", \"$OUTPUT\")"
done

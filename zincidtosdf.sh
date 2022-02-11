i=$(grep -n $1 16/zincid.txt | cut -d: -f1)
is=$((i-2))
ie=$((i-1))
cs=$(tail -c +$((1+8*is)) 16/mconfs.u64 | head -c 8 | decodeu64)
ce=$(tail -c +$((1+8*ie)) 16/mconfs.u64 | head -c 8 | decodeu64)
os=$(tail -c +$((-7+8*cs)) 16/ligand.ftr | head -c 8 | decodeu64)
oe=$(tail -c +$((-7+8*ce)) 16/ligand.ftr | head -c 8 | decodeu64)
tail -c +$((1+os)) 16/ligand.sdf | head -c $((oe-os)) > $1.sdf

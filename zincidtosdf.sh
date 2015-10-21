i=$(grep -n $1 16_zincid.txt | cut -d: -f1)
is=$((i-2))
ie=$((i-1))
cs=$(tail -c +$((1+8*is)) 16_mconfs.u64 | head -c 8 | decodeu64)
ce=$(tail -c +$((1+8*ie)) 16_mconfs.u64 | head -c 8 | decodeu64)
os=$(tail -c +$((-7+8*cs)) 16_ligand.ftr | head -c 8 | decodeu64)
oe=$(tail -c +$((-7+8*ce)) 16_ligand.ftr | head -c 8 | decodeu64)
tail -c +$((1+os)) 16_ligand.sdf | head -c $((oe-os)) > $1.sdf

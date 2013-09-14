#rm -f snapshot_T*.txt*
#Ising_SA 50 0.9999 1000
foreach file (`ls snapshot_T*.txt`)
plot_snapshot snap.cpt $file region 1
end

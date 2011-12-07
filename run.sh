NODE_MEMORY=`free -k | grep ^Mem: | awk '{ print $2; }'`
NODE_MEMORY_LIMIT=`echo "0.95 * $NODE_MEMORY / 1" | bc`
#ulimit -v $NODE_MEMORY_LIMIT -m $NODE_MEMORY_LIMIT
echo "memory limit: $NODE_MEMORY_LIMIT kilobytes"
#ibrun ./reference_glm -f ./reference_glm.massive -v 2
#-g /home/01722/tpatki/iplant/fstest/fstest/massiver.dat \
#-g /scratch/01713/pbailey/massiver.dat \
#-f /scratch/01713/pbailey/data_maize_nam/fixed.effects.nam.sorted.filtered.bin \
ibrun ./reference_glm \
-f /scratch/01722/tpatki/fixed.effects.nam.sorted.filtered.bin \
--num_fixed 26 \
-g /scratch/01722/tpatki/single.dat \
--num_geno 500000 \
-r /scratch/01722/tpatki/residuals.chr10.sorted.bin \
--num_r 4892 \
-v2 \
$1

NODE_MEMORY=`free -k | grep ^Mem: | awk '{ print $2; }'`
NODE_MEMORY_LIMIT=`echo "0.95 * $NODE_MEMORY / 1" | bc`
#ulimit -v $NODE_MEMORY_LIMIT -m $NODE_MEMORY_LIMIT
echo "memory limit: $NODE_MEMORY_LIMIT kilobytes"
#ibrun ./reference_glm -f ./reference_glm.massive -v 2

ibrun ../reference_glm \
-f /scratch/01713/pbailey/data_maize_nam/fixed.effects.nam.sorted.filtered.bin \
--num_fixed 26 \
-g /scratch/01713/pbailey/4892x1M_double.dat \
--num_geno 1000000 \
-r /scratch/01713/pbailey/data_maize_nam/residuals.chr10.sorted.bin \
--num_r 4892 \
-v2 \
$1 > log

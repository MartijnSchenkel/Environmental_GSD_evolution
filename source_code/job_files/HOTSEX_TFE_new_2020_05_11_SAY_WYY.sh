#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-00:00:00
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000MB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.a.schenkel@rug.nl

# prep work
DATE=2020_05_11_SAY_WYY
HSDIR=$HOME/2019_02_HOTSEX/source_min_parameters
module load GCC/8.2.0-2.31.1

mkdir $HOME/2019_02_HOTSEX/data/$DATE
mkdir /data/p275703/$DATE

for n in {8000..9999}
do

SIMDIR=$HOME/2019_02_HOTSEX/data/$DATE/$n
g++ -I $SIMDIR -x c++ $HSDIR/*.cpp -x c++-header $HSDIR/*.h -x c++-header $SIMDIR/parameters.h -o $SIMDIR/HOTSEX -std=c++17 -O2

cat << EOT >> TFE$DATE$n.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=3-00:00:00
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000MB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.a.schenkel@rug.nl

cd $SIMDIR
$SIMDIR/HOTSEX

cd /home/p275703
mv $SIMDIR /data/p275703/$DATE
EOT

sbatch TFE$DATE$n.sh
rm TFE$DATE$n.sh

done

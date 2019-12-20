#!/bin/bash 

#SBATCH --nodelist=europa

MIN=0.0
MAX=2.0
STP=0.1

rm -vf tra.txt
rm -vf rdw.txt
rm -vf lam.txt

sed -i "s/lagrange_multiplier = xx/lagrange_multiplier = 1/g" etc/phynet.in
sed -i "s/loss_type = lt/loss_type = rd/g" etc/phynet.in
for rdw in $(seq $MIN $STP $MAX)
do
	sed -i "s/random_domain_bound = xx/random_domain_bound = $rdw/g" etc/phynet.in
   	phynet-silent.x | tail -1 >> rdw.txt
   	sed -i "s/random_domain_bound = $rdw/random_domain_bound = xx/g" etc/phynet.in
   	tail -1  rdw.txt
done
sed -i "s/loss_type = rd/loss_type = lt/g" etc/phynet.in

sed -i "s/loss_type = lt/loss_type = pg/g" etc/phynet.in
for tra in $(seq $MIN $STP $MAX)
do
	sed -i "s/trade_off_parameter = xx/trade_off_parameter = $tra/g" etc/phynet.in
   	phynet-silent.x | tail -1 >> tra.txt
   	sed -i "s/trade_off_parameter = $tra/trade_off_parameter = xx/g" etc/phynet.in
   	tail -1  tra.txt
done
sed -i "s/loss_type = pg/loss_type = lt/g" etc/phynet.in


sed -i "s/lagrange_multiplier = 1/lagrange_multiplier = xx/g" etc/phynet.in
sed -i "s/loss_type = lt/loss_type = c2/g" etc/phynet.in
for lambda in $(seq $MIN $STP $MAX)
do
	sed -i "s/lagrange_multiplier = xx/lagrange_multiplier = $lambda/g" etc/phynet.in
	phynet-silent.x | tail -1 >> lam.txt
	sed -i "s/lagrange_multiplier = $lambda/lagrange_multiplier = xx/g" etc/phynet.in
	tail -1  lam.txt
done
sed -i "s/loss_type = c2/loss_type = lt/g" etc/phynet.in


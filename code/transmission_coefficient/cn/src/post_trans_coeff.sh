echo 'import numpy as np
import sys


def trans_coeff(filename):
    table  = np.loadtxt(filename)
    # initial step position
    pos0   = table[0,:]
    # initial step velocity
    vel0   = table[1,:] - pos0
    #print(vel0)
    r_hat_table = table - pos0
    r_hat_table = r_hat_table[1:,:]
    # step function theta:
    # for positive direction: theta = 1 if r_hat > 0
    #                         theta = 0 if r_hat < 0
    theta0 = np.zeros_like(vel0)
    theta0[np.sign(vel0) == -1] = 1
    theta_table = np.zeros_like(r_hat_table)
    theta_table[ np.sign(r_hat_table) == -1 ] = 1
    #print(theta_table[0:2])
    # initial step flux, averaged
    j0_pos_avg = np.sum(vel0 * theta0     ) / len(vel0)
    # flux at time t, averaged
    jt_pos_avg = np.sum(vel0 * theta_table, axis = 1) / len(vel0)
    kappa_t = jt_pos_avg / j0_pos_avg
    #print(j0_pos_avg)
    #print(jt_pos_avg)
    time_step = range(len(kappa_t))
    [print(t,kt) for t,kt in zip(time_step,kappa_t)]



if __name__ == "__main__":
    filename=str(sys.argv[1])
    trans_coeff(filename)
    exit()

' > trans_coeff.py
mkdir -p post_process
echo "extra CV time series from individual folders" 
echo "coord num is the third column"
k=0 
for i in $(ls -d config_*)
do 
	head -n 1002 $i/COLVAR-comm |grep -v '#' |awk '{print $3}' > post_process/$i.dat &
	k=$(($k+1))
	if [[ $k -eq 10 ]]; then wait; fi
done
wait
echo "combining CV time series" && paste post_process/config_*.dat > post_process/cv_table

#echo "remove temp files" && rm post_process/config_* &

echo "running python calculations" &&\

timestep=0.001
python trans_coeff.py post_process/cv_table | awk -v dt=${timestep} '{printf("%14.10f\t%14.10f\n",$1*dt,$2)}' > trans_coeff.dat

wait

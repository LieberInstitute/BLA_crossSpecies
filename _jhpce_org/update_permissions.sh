#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=400G
#$ -N update_permissions
#$ -o logs/update_permissions.txt
#$ -e logs/update_permissions.txt


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module list

MAINDIR="/dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/"

## Remove default permissions
find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -x "A:g:lieber_lcolladotor@cm.cluster:rwaDxtcy" {} \;
find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -x "A:gfdi:lieber_lcolladotor@cm.cluster:rwaDxtcy" {} \;
find ${MAINDIR} -user ${USER} -type f -exec nfs4_setfacl -x "A:g:lieber_lcolladotor@cm.cluster:rwatcy" {} \;

find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -x "A:g:lieber_marmaypag@cm.cluster:rwaDxtcy" {} \;
find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -x "A:gfdi:lieber_marmaypag@cm.cluster:rwaDxtcy" {} \;
find ${MAINDIR} -user ${USER} -type f -exec nfs4_setfacl -x "A:g:lieber_marmaypag@cm.cluster:rwatcy" {} \;

## Set new permissions
find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RWX" {} \;
find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -a "A:gfdi:hickslab@cm.cluster:RWX" {} \;
find ${MAINDIR} -user ${USER} -type f -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RW" {} \;

find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RWX" {} \;
find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_lcolladotor@cm.cluster:RWX" {} \;
find ${MAINDIR} -user ${USER} -type f -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RW" {} \;

find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RWX" {} \;
find ${MAINDIR} -user ${USER} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_marmaypag@cm.cluster:RWX" {} \;
find ${MAINDIR} -user ${USER} -type f -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RW" {} \;


## Make sure files are under the right user group
chgrp lieber_marmaypag -R ${MAINDIR}

## And the right permissions
chmod -R 774 ${MAINDIR}

## For setting the group sticky bit
find ${MAINDIR} -type d | xargs chmod g+s

nfs4_getfacl ${MAINDIR}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

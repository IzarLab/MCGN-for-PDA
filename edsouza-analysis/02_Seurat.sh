#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda deactivate
conda activate cxcr4

for i in {18..40}
do
    if [[  $i -ne 25 ]]
    then
        pat="GM${i}"

        outdir="/home/ubuntu/edsouza-summer2023/cxcr4-pdac/data/seurat/${pat}"
        logfile="${outdir}/log.txt"
        donefile="${outdir}/finished.checkpoint"

        if [ ! -f "${donefile}" ]
        then
            echo "Starting ${pat}" | ts '[%Y-%m-%d %H:%M:%S]'
            mkdir -p "${outdir}"
            Rscript /home/ubuntu/edsouza-summer2023/cxcr4-pdac/02b_Seurat.R $i \
                2>&1 | ts '[%Y-%m-%d %H:%M:%S]' \
                    > "${logfile}" 
                    
            if ! grep -q "Execution halted" <(tail -n 2 "${logfile}"); 
            then
                echo "Sample ${pat} finished successfully"
                touch "${donefile}"
            else
                echo "Sample ${pat} finished unsuccessfully"
            fi

            echo "Finished ${samplenum}" | ts '[%Y-%m-%d %H:%M:%S]'
        fi

    fi
done
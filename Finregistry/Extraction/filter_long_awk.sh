zcat /home/mattferr/Projects/SD-Connect/project_2007099/processed_data/FinRegistry_v01/detailed_longitudinal/detailed_longitudinal_2024-24-01.csv.gz | awk -v OFS=',' '
BEGIN { FS = OFS; }
NR == FNR { 
    # Load IDs from the ID file into an array
    id[$1] = 1; 
    next; 
}
# Now process the main CSV data
$1 in id && $2 ~ /^(INPAT|OUTPAT|OPER_IN|OPER_OUT|PRIM_OUT|CANC|PURCH)$/ {
    print $1, $2, $3, $4, $5, $14;
}
' /media/volume/mferro/scripts/CCI/idfile_adherence_all.csv - | gzip > /media/volume/mferro/data/CCI/longitudinal_adherence_statins_all.csv.gz

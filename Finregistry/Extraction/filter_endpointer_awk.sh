zcat /home/mattferr/Projects/SD-Connect/project_2007099/processed_data/FinRegistry_v00/endpointer/densified_first_events_DF10_no_omits_2022-09-20.csv.gz | awk -v OFS=',' '
BEGIN { FS = OFS; }
NR == FNR { 
    # Load IDs from the ID file into an array
    id[$1] = 1; 
    next; 
}
# Now process the main CSV data
$1 in id  {
    print $1, $2, $3;
}
' /media/volume/mferro/scripts/CCI/idfile_adherence_all.csv - > /media/volume/mferro/data/STATIN_USERS_ENDPOINTS.csv

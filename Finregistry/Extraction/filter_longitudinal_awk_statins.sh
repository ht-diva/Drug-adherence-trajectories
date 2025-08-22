

gzip -dc /home/mattferr/Projects/SD-Connect/project_2007099/processed_data/FinRegistry_v01/detailed_longitudinal/detailed_longitudinal_2024-24-01.csv.gz | awk 'BEGIN { 
	FS = OFS = ","  # Set field separator to comma
	cmd = "gzip -c > /media/volume/mferro/data/DRUGS/STATINS.csv.gz" 
}

  # For the large CSV input file (inputfile.csv)
  ($5 ~ /^C10AA/ || $5 ~ /^C10B/) {print | cmd
  } END  {
	close(cmd) }' 




import gzip

WORKING_DIR = "/stanley/huang_lab/home/myu/IBD/preprocess"
manifest_file = WORKING_DIR + "/pheno_data/gnomADv4.manifest.03-02.txt" # gs://ibd-exomes-gnomad-subset/rzhu/00.data/gnomADv4.manifest.03-02.txt

IID2imputed_sex = {}

line_count = 0
for line in gzip.open(WORKING_DIR + "/QC/sex_check/gnomad_ibd_imputed_sex.tsv.bgz", "rt"): # upload to gs://ibd-exomes-gnomad-subset/QC_round2/6.export/gnomad_ibd_imputed_sex.tsv.bgz
    parsed_line = line.strip('\n').split("\t")
    if line_count == 0:
        col2idx = {}
        for i in range(len(parsed_line)):
            col2idx[parsed_line[i]] = i
    else:
        IID = parsed_line[col2idx['s']].replace(" ", ",")
        if parsed_line[col2idx['is_female']] == "true":
            IID2imputed_sex[IID] = "Female"
        elif parsed_line[col2idx['is_female']] == "false":
            IID2imputed_sex[IID] = "Male"
        else:
            IID2imputed_sex[IID] = "NA"
            
    line_count += 1
    
line_count = 0
for line in gzip.open(WORKING_DIR + "/QC/sex_check/moayeddi_imputed_sex.tsv.bgz", "rt"):
    parsed_line = line.strip('\n').split("\t")
    if line_count == 0:
        col2idx = {}
        for i in range(len(parsed_line)):
            col2idx[parsed_line[i]] = i
    else:
        IID = "RP-1915_" + parsed_line[col2idx['s']].replace(" ", ",").replace("\"", "")
        if parsed_line[col2idx['is_female']] == "true":
            IID2imputed_sex[IID] = "Female"
        elif parsed_line[col2idx['is_female']] == "false":
            IID2imputed_sex[IID] = "Male"
        else:
            IID2imputed_sex[IID] = "NA"
    line_count += 1

output = open(WORKING_DIR + "/QC/sex_check/all.sex_check.txt", "w+")
output.write("\t".join(['IID', 'reported_sex', 'imputed_sex', 'sex_final']) + "\n")
line_count = 0
for line in open(manifest_file, "r"):
    parsed_line = line.strip('\n').split("\t")
    if line_count == 0:
        col2idx = {}
        for i in range(len(parsed_line)):
            col2idx[parsed_line[i]] = i
    else:
        IID = parsed_line[col2idx['analysis_id']].replace("\"", "").replace(",MS", "")
        reported_sex = parsed_line[col2idx['sex']] if parsed_line[col2idx['sex']] in ['Male', 'Female'] else "NA"
        imputed_sex = IID2imputed_sex[IID] if IID2imputed_sex.get(IID) else "NA"
        
        if reported_sex in ['Male', 'Female'] and imputed_sex in ['Male', 'Female']:
            if parsed_line[col2idx['sex']] == IID2imputed_sex[IID]:
                sex_final = reported_sex
            else:
                sex_final = "Inconsistent"
        elif not reported_sex in ['Male', 'Female'] and not imputed_sex in ['Male', 'Female']:
            sex_final = "NA"
        elif reported_sex in ['Male', 'Female']:
            sex_final = reported_sex
        else:
            sex_final = imputed_sex
        
        output.write("\t".join([IID, reported_sex, imputed_sex, sex_final]) + "\n")
        
    line_count += 1

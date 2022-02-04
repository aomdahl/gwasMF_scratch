import glob
joined = set()
files = glob.glob("/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/*allpairs.txt")
for file in files:
    print(file)
    with open(file, 'r') as istream:
        counter = 0
        for l in istream:
            counter += 1
            line = l.strip().split()[1]
            #chr1_261402_CTT_C_b38
            snp = line[0:-4].replace("_", ":").replace("chr", "")
            joined.add(snp)
            if(counter > 1000000):
                print("one mil in")
                counter = 0
with open("/work-zfs/abattle4/ashton/reference_data/GTEx_v8_snps.all.txt", 'w') as ostream:
    for snp in joined:
        ostream.write(snp + "\n")

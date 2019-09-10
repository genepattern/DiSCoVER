# docker run --rm -w /job_1234 -v $PWD/job_1234:/job_1234 -v $PWD/data:/temp/data -it genepattern/discover:1.0 ls /module

# to test with double stranded data
# docker run --rm -w /job_1234 -v $PWD/job_1234:/job_1234 -v $PWD/data:/temp/data -it genepattern/discover:1.0 python /module/call_discover.py -g "/temp/data/input_file_list.txt" -m True -u False -c "nothing to be seen here"

# docker run --rm -w /job_1234 -v $PWD/job_1234:/job_1234 -v $PWD/data:/temp/data -it genepattern/discover:1.0 python /module/call_discover.py -g "https://datasets.genepattern.org/data/module_support_files/DiSCoVER/gene_abundance.sleuth.csv" -m True -u False -c "nothing to be seen here"

docker run --rm -w /job_1234 -v $PWD/job_1234:/job_1234 -v $PWD/data:/temp/data -it genepattern/discover:1.3 python /module/call_discover.py -g "/temp/data/gene_expression.csv" -u "Cerebellar" -c "nothing to be seen here"

# docker run --rm -w /job_1234 -v $PWD/job_1234:/job_1234 -v $PWD/data:/temp/data -it genepattern/discover:1.1 python /module/call_discover.py -g "https://datasets.genepattern.org/data/module_support_files/DiSCoVER/gene_abundance.sleuth.csv" -m True -u False -c "nothing to be seen here"

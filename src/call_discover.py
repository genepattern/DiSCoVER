#==============================================================================
#NOTE: If you are looking for the code of DiSCoVER itself, read the README and then head to: https://datasets.genepattern.org/?prefix=data/module_support_files/DiSCoVER/
#==============================================================================
print('About to start running the module')
### GP Module reqs
import argparse
from timeit import default_timer as timer
import shutil
beginning_of_time = timer()

parser = argparse.ArgumentParser()
# ~~~~Module Required Arguments~~~~~ #
parser.add_argument("-g", "--gene_expression",
                    type=str,
                    help="Name of the csv file which contains the gene expression to be used")
# parser.add_argument("-m", "--medullo",
#                     type=str,
#                     help="Whether or not this sample classified as medulloblastoma",
#                     default='False')
parser.add_argument("-u", "--use_control",
                    type=str,
                    help="What type of control to use",
                    default='Custom')
parser.add_argument("-c", "--control",
                    type=str,
                    help="Name of the csv file which contains the gene expression of the control",
                    default='')
parser.add_argument("-s", "--supplementary",
                    type=str,
                    help="Whether or not to output supplementary files",
                    default='False')

# ~~~~Development Optional Arguments~~~~~ #
# Reminder: "store_true" args, are False by default and when the flag is passed
# then they become True
parser.add_argument("-v", "--verbose",
                    action="store_true",
                    help="increase output verbosity")
parser.add_argument("-d", "--debug",
                    action="store_true",
                    help="increase output verbosity")
args = parser.parse_args()
if args.verbose:
    print("Ah! The old verbosaroo")

print("~~~~~~~~~~~~~~~~~~~~~~")
print("Using arguments:")
print(args)
print("Now getting work done.")
print("~~~~~~~~~~~~~~~~~~~~~~")

def read_list_of_files(filename):
    counter = 0
    patient_exp_file = ''
    with open(filename) as file:
        for line in file.readlines():
            patient_exp_file = line.strip('\n')
            counter +=1
    if counter > 1:
        log('WARNING: There is more than one file provided, this is unnexpected!')
    return patient_exp_file

### end of GP module requirements

# from discover import *
from companion_script import *

# from rpy2.robjects import numpy2ri
# numpy2ri.activate()
from discover import discover_from_expression, plot_discover_from_expression
from drug_suggestion.expression.controls import load_control_exp

# patient_exp = pd.read_csv('test_data/gene_abundance.sleuth.csv', index_col=0)
# patient_exp = pd.read_csv('RNASeq_quant/gene_expression.csv', index_col=0)
# patient_exp = pd.read_csv('RNASeq_quant/gene_abundance.sleuth.csv', index_col=0)

# This assumes only one file is provided
patient_exp = pd.read_csv(args.gene_expression,index_col=0)

if len(patient_exp)>1:
    # This means that the csv file was tall, not wide (as expected):
    patient_exp = patient_exp.set_index(patient_exp.columns[0]).T

#TODO: Set this as a parameter
# patient_exp.index.values[0]
case_id = patient_exp.index.values[0]

# Parse control
if args.use_control=='Custom':
    expression_control = 'custom_control'
    control_exp = pd.read_csv(args.control,index_col=0)
    if len(control_exp)>1:
        # This means that the csv file was tall, not wide (as expected):
        control_exp = control_exp.set_index(control_exp.columns[0]).T
    log(f'Using the custom control provided by the file {read_list_of_files(args.control)}')
elif args.use_control=='Cerebellar':
    expression_control = 'cerebellar_stem'
    control_exp = load_control_exp(f'/build/drug_suggestion/expression/controls/{expression_control}')
    log(f'Using the built-in control {expression_control}.')
elif args.use_control=='Neural':
    expression_control = 'neural_stem'
    control_exp = load_control_exp(f'/build/drug_suggestion/expression/controls/{expression_control}')
    log(f'Using the built-in control {expression_control}.')
else:
    exit(f'Unexpected value for parameter "use_control", value {args.use_control}')

discover_out_dir = 'supplementary_files'
full_discover_results_file = os.path.join(discover_out_dir,'discover.all.csv')

log("About to perform DiSCoVER.")

discover_results = discover_from_expression(exp=patient_exp,
                                            control_exp=control_exp,
                                            verbose=False, extra_outputs=True)
log("DiSCoVER has finished. Now onto some post processing")

# move some files created by DiSCoVER
if not os.path.exists(discover_out_dir):
    os.mkdir(discover_out_dir)

for current_file in ['cell_lines_IDs_and_types_ccle.csv','cell_lines_IDs_and_types_COSMIC_IDS_gdsc.csv','cell_lines_IDs_and_types_ctrp.csv']:
        os.rename(current_file,os.path.join(discover_out_dir,current_file))
        log(f'Moved {current_file} to {os.path.join(discover_out_dir,current_file)}')

log("Ranking cell lines by enrichment and saving those")

temp_setup = {}
temp_setup['case_id'] = case_id
temp_setup['discover_out_dir'] = discover_out_dir
setup = Bunch(temp_setup)

ranked_diseases_from_enrichment = rank_cell_lines(setup)
log('Saving results to file.')
discover_results.T.sort_values(by=setup.case_id, ascending=False).to_csv(full_discover_results_file)
log("Saving done!")


from drug_suggestion.drug_annotation import subset_to_reasonable_drugs
from drug_suggestion.expression.discover import load_discover_drug_to_cids
disco2cid = load_discover_drug_to_cids()
reasonable_results = subset_to_reasonable_drugs(discover_results,
                                            disco2cid,
                                            out_prefix='discover.{}'.format(expression_control),
                                            out_dir=discover_out_dir)

# This will override the file setup.rdrugs_discover_file
all_drugs = format_drugs(discover_results,
                            disco2cid,
                            out_prefix='discover.{}'.format(expression_control),
                            out_dir=discover_out_dir)

if args.supplementary=='False':
    shutil.rmtree(discover_out_dir)

format_drugs(discover_results,
             disco2cid,
             out_prefix='discover.{}'.format(expression_control),
             out_dir='.')

log("DiSCoVER done!")


end_of_time = timer()
print("We are done! Wall time elapsed (in minutes):", (end_of_time - beginning_of_time)/60)

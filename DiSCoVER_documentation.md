# DiSCoVER (v1.0)

This is a bare-bones implementation of DiSCoVER, intended to be used for drug recommendation based on RNAseq expression (e.g., the results from Kallisto).

Author: Edwin Juarez

Contact: https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help

Algorithm Version: DiSCoVER 1.0

## Summary
*To be added*

## References

DiSCoVER paper:
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055054/
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055054/?report=reader
- http://clincancerres.aacrjournals.org/content/22/15/3903
- Hanaford, A. R., Archer, T. C., Price, A., Kahlert, U. D., Maciaczyk, J., Nikkhah, G., … Raabe, E. H. (2016). DiSCoVERing Innovative Therapies for Rare Tumors: Combining Genetically Accurate Disease Models with In Silico Analysis to Identify Novel Therapeutic Targets. Clinical cancer research : an official journal of the American Association for Cancer Research, 22(15), 3903–3914. doi:10.1158/1078-0432.CCR-15-3011

*More to be added*

### Functionality yet to be implemented:
*To be added*

### Technical notes:
*To be added*

## Parameters

### -->gene_expression
type=str  
The csv file which contains the gene expression to be used  

### -->is_medulloblastoma
type=str (choice: 'True' or 'False')  
help="Whether or not this sample classified as medulloblastoma"  
default='False'

### -->use_custom_control
type=str (choice: 'True' or 'False')
Whether or not to use a custom control  
default='False'  
### -->control
type=str,  
The csv file which contains the gene expression of the control"  
default=''  

## License

DiSCoVER itself and this GenePattern module are distributed under a modified BSD license. This module's license is available at https://github.com/genepattern/DiSCoVER/blob/master/LICENSE

## Platform Dependencies
Task Type: Drug recommendation
CPU Type:
any

Operating System:
any

Language:
python 3.7

Version Comments
Version	Release Date	Description
1	2019-05-21	Initial release of DiSCoVER

# Spy-C predictions

1. File contains the following columns
- gene
- Uniprot
- PTM site
- oriented_peptide
- 6mer_peptide : Input for SpY-C
- predicted_class : Binder (1); Nonbinder (0)
- binder_prob : predicted probability
- confidence_category : >=0.6 is confident binder; <=0.4 is confident nonbinder; anything else is 'low-confidence'


2. If a peptide has nan values for predicted_class, binder_prob and confidence_category -> short peptide (ones towards the N/C terminii)
3. If a peptide has nan values for predicted_class and binder_prob -> is found in SpY-C training dataset (their ground truth label in 'confidence_category' column)

   


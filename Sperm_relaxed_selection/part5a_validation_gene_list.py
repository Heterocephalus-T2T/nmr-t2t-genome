#!/usr/bin/env python3
"""Part 5a: write the validation gene list and its scheme config.

The 17 control genes (8 positive, 9 negative) with their literature-based
test species, curated in the table below. Two files come out:

  gene_lists/validation.tsv   SYMBOL, TYPE, EXPECTED_K, TEST_SPECIES,
                              SPECIFIC_REFERENCE, CITATION   (the record)
  schemes/validation.tsv      Scheme, Gene, Test_Species, Reference_Species
                              (the part-9 contract; general mode kept per the
                              report: Reference = ALL_OTHERS)

Species use pipeline codes (NMR = HetGla, BATS = BatSui, AMR = FukAns).

Usage: python3 part5a_validation_gene_list.py
"""

import os

# gene: (type, expected_K, test_species, specific_reference, citation)
CONTROLS = {
    'CRYAA':   ('positive', '<1',  'NMR,BATS,AMR', 'MOUS,RAT,HUMA,MACA,DOG,CAT',
                'Emerling & Springer 2014'),
    'CRYAB':   ('positive', '<1',  'NMR,BATS,AMR', 'MOUS,RAT,HUMA,MACA,DOG,CAT',
                'Emerling & Springer 2014'),
    'CRYBB1':  ('positive', '<1',  'NMR,BATS,AMR', 'MOUS,RAT,HUMA,MACA',
                'Emerling & Springer 2014'),
    'RHO':     ('positive', '<=1', 'NMR,BATS,AMR', 'MOUS,RAT,HUMA,MACA,DOG,CAT',
                'Emerling & Springer 2014'),
    'GRK1':    ('positive', '<1',  'NMR,BATS,AMR', 'MOUS,RAT,HUMA,MACA',
                'Emerling & Springer 2014'),
    'OPN1SW':  ('positive', '<1',  'NMR,BATS',     'SQUI,HUMA,MACA,GPIG,CHIN,BONO,GORI',
                'Zhao et al. 2009'),
    'OR1A1':   ('positive', '<1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,DOG,CAT,ELEP,PIG,BOS',
                'Gilad et al. 2004'),
    'OR2W1':   ('positive', '<1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,DOG,CAT,ELEP,PIG,BOS',
                'Gilad et al. 2004'),
    'ATP5F1A': ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'de Jonge 2007; Eisenberg & Levanon 2013'),
    'RPL13':   ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'de Jonge 2007; Eisenberg & Levanon 2013'),
    'RPS13':   ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'de Jonge et al. 2007'),
    'OPN1LW':  ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'Zhao et al. 2009; Jacobs 2013'),
    'VCP':     ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'Eisenberg & Levanon 2013'),
    'GPI':     ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'Eisenberg & Levanon 2013'),
    'REEP5':   ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'Eisenberg & Levanon 2013'),
    'VPS29':   ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'Eisenberg & Levanon 2013'),
    'SNRPD3':  ('negative', '~1',  'HUMA,MACA,BONO,GORI', 'MOUS,RAT,GPIG,CHIN,DEGU',
                'Eisenberg & Levanon 2013'),
}


def main():
    os.makedirs('gene_lists', exist_ok=True)
    os.makedirs('schemes', exist_ok=True)

    with open('gene_lists/validation.tsv', 'w') as f:
        f.write("SYMBOL\tTYPE\tEXPECTED_K\tTEST_SPECIES\tSPECIFIC_REFERENCE\tCITATION\n")
        for g, (typ, ek, test, ref, cit) in CONTROLS.items():
            f.write(f"{g}\t{typ}\t{ek}\t{test}\t{ref}\t{cit}\n")

    with open('schemes/validation.tsv', 'w') as f:
        f.write("Scheme\tGene\tTest_Species\tReference_Species\n")
        for g, (typ, ek, test, ref, cit) in CONTROLS.items():
            f.write(f"general\t{g}\t{test}\tALL_OTHERS\n")

    n_pos = sum(1 for v in CONTROLS.values() if v[0] == 'positive')
    n_neg = len(CONTROLS) - n_pos
    print(f"validation set: {len(CONTROLS)} genes ({n_pos} positive, {n_neg} negative)")
    print("wrote gene_lists/validation.tsv and schemes/validation.tsv")


if __name__ == '__main__':
    main()

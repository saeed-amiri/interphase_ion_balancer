"""
Reading .gro file
Goal: push ions from interface to an area below the nanoparticle (NP)
How:
    1- get gro file
    2- find nanoaprticle depth
    3- define the place of the wall based on the NP depth
    4- find the index of the ions above that wall
    5- find posioions for that number of ions in the area below NP
    6- update the positions of the selected ions with the new positions
    7- needs an input file, to set the input values:
        - Names of the ions (list)
        - Names of the NPs (list)
        - Other inputs from update_strucure
Note:
    1- Keep in mind the douple NP system!
    2- The condition for putting several wall!
"""

import sys
import logger
import update_residues_gro
import write_gro_file


if __name__ == "__main__":
    GRO_ATOMS = update_residues_gro.UpdateResidues(
        sys.argv[1], log=logger.setup_logger('balancer.log'))
    write_gro_file.write_gromacs_gro(gro_data=GRO_ATOMS.combine_residues,
                                     filename='system_ion_balanced.gro',
                                     title=GRO_ATOMS.title,
                                     pbc_box=GRO_ATOMS.pbc_box)

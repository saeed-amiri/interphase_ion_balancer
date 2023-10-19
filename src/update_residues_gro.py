"""To properly format GRO files, velocities columns must be included
in the input file. 
"""

import sys
import pandas as pd
import logger
import ionization
from colors_text import TextColor as bcolors


class UpdateResidues:
    """get all the dataframes as an object
        and update the position of the selected ions
    """

    # Message to pass for logging and writing
    info_msg: str = 'Messages from UpdateResidues:\n'
    # To append all the residues to the dict
    updated_residues: dict[str, pd.DataFrame] = {}
    # To append new atoms to the main atoms
    updated_atoms: pd.DataFrame
    # To return new atoms seperatly
    new_ions: pd.DataFrame
    # Final datafrme
    combine_residues: pd.DataFrame
    nr_atoms_residues: dict[str, dict[str, int]] = {}
    # Title and pbc of the system to write in the final file
    title: str
    pbc_box: str

    def __init__(self,
                 fname: str,  # Name of the input file (pdb)
                 log
                 ) -> None:
        data = ionization.IonizationSol(fname, log)
        self.title = data.title
        self.pbc_box = data.pbc_box
        self.param = data.param
        data.residues_atoms['CLA'] = \
            self.update_ions(data.residues_atoms['CLA'].copy(),
                             data.updated_ions.copy())
        combine_residues = self.concate_residues(data.residues_atoms)
        self.combine_residues = self.__set_atom_id(combine_residues,
                                                   data.param['DEBUG'])
        self.write_log_msg(log)

    def update_ions(self,
                    all_ions: pd.DataFrame,
                    updated_ions: pd.DataFrame
                     ) -> None:
        """replace the up_wall_ions with updated ones"""
        # Set 'residue_number' as the index in both DataFrames
        for residue in updated_ions.index:
            all_ions.loc[residue] = updated_ions.loc[residue]
        return all_ions

    def concate_residues(self,
                         updated_residues: dict[str, pd.DataFrame]
                         ) -> pd.DataFrame:
        """concate all the residues in one dataframe, The order is
        very important. it should follow the order in the main file"""
        cols_order: list[str] = \
            ['SOL', 'D10', 'ODN', 'CLA', 'POT', 'ODM', 'COR', 'APT']
        # Drop unesist residues
        for res in ['ODN', 'POT', 'ODM']:
            if res not in updated_residues:
                cols_order.remove(res)
        # Concatenate DataFrames in the desired order
        combine_residues: pd.DataFrame = \
            pd.concat([updated_residues[col] for col in cols_order],
                      axis=0, ignore_index=True)

        return combine_residues

    @staticmethod
    def __set_atom_id(combine_residues: pd.DataFrame,  # All the rsidues
                      debug: bool
                      ) -> pd.DataFrame:
        """set the atom_id for the all the atoms"""
        df_c: pd.DataFrame = combine_residues.copy()
        # Specify the limit for the atom IDs
        atom_id: list[int] = \
            mk_atom_id_cycle(len(combine_residues), start_id=1)
        # Calculate the number of cycles
        df_c['atom_id'] = atom_id
        if debug != 'None':
            df_c.to_csv('combine_residues.debug', sep=' ')
        return df_c

    @staticmethod
    def get_atoms(atoms: pd.DataFrame,  # Initial system
                  new_hn3: pd.DataFrame,  # New NH3 atoms
                  new_ions: pd.DataFrame,  # New Ions atoms
                  ) -> pd.DataFrame:
        """append the new atoms to the main dataframe with all atoms"""
        return pd.concat([atoms, new_hn3, new_ions])

    def write_log_msg(self,
                      log: logger.logging.Logger  # Name of the output file
                      ) -> None:
        """writing and logging messages from methods"""
        log.info(self.info_msg)
        print(f'{bcolors.OKBLUE}{UpdateResidues.__module__}:\n'
              f'\t{self.info_msg}\n{bcolors.ENDC}')


# Helper function to update index in gro fasion
def mk_atom_id_cycle(list_len: int,  # Size of the list,
                     start_id: int,  # First index of the residues
                     id_limit=99999  # Limit of the cycle
                     ) -> list[int]:
    """
    Generate a list of unique atom IDs in a custom cycle.

    This function generates a list of unique integers that follows a
    custom cycle pattern for atom IDs. The first cycle starts from the
    provided 'start_id' and goes up to 'id_limit', and subsequent
    cycles continue from 0 to 'id_limit'.

    Parameters:
        list_len (int): The desired size of the list to be generated.
                        The function will generate unique atom IDs
                        until reaching this size.
        start_id (int, optional): The starting value for the first
                                  cycle. The default value is 1.
        id_limit (int, optional): The upper limit of the atom ID cycle.
                                  The default value is 99999.

    Returns:
        List[int]: A list of unique atom IDs in the custom cycle.

    Example:
        >>> mk_atom_id_cycle(5)
        [1, 2, 3, 4, 5]
        >>> mk_atom_id_cycle(8, start_id=10)
        [10, 11, 12, 13, 14, 0, 1, 2]
        >>> mk_atom_id_cycle(12, start_id=100, id_limit=105)
        [100, 101, 102, 103, 104, 105, 0, 1, 2, 3, 4, 5]
    """
    counter = 0
    atoms_id = []
    while counter < list_len:
        if counter == 0:
            cycle_i = start_id
            cycle_f = id_limit + 1 - start_id
        else:
            cycle_i = 0
            cycle_f = id_limit + 1
        slice_i = [item + cycle_i for item in range(cycle_f)]
        counter += len(slice_i)
        atoms_id.extend(slice_i)
        del slice_i
    return atoms_id[:list_len]


if __name__ == '__main__':
    UpdateResidues(sys.argv[1], log=logger.setup_logger('update.log'))

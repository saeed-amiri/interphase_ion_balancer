"""
ProcessData Class
-----------------

The ProcessData class is designed to efficiently process data from a
pdb or gro file and extract relevant sections of data for different
residues or residue groups. It creates separate data frames for diffe-
rent residues or groups of residues and performs various operations to
analyze and manipulate the data. Calculating the diameter of nanopart-
icles (NPs) basedv on APTES positions.

Attributes:
    atoms (pd.DataFrame): A pandas DataFrame containing data for all
            atoms in the system.
    param (dict[str, float]): A dictionary containing various paramet-
            ers from an input file.
    residues_atoms (dict[str, pd.DataFrame]): A dictionary containing
            pandas DataFrames for each residue's atoms.
    np_diameter (np.float64): The maximum radius of the nanoparticle
            (NP) based on APTES positions.
    title (str): The name of the system; applicable if the file is in
            gro format.
    pbc_box (str): The periodic boundary condition of the system;
            applicable if the file is in gro format.

Methods:
    __init__(fname: str, log: logger.logging.Logger) -> None:
        Initialize the ProcessData object.

    calculate_maximum_np_radius() -> np.float64:
        Calculate the maximum radius of the nanoparticle (NP) based on
        APTES positions.

    get_unique_residue_names() -> list[str]:
        Get the list of unique residue names in the system.

Private Methods:
    _get_data(fname: str, log: logger.logging.Logger) -> pd.DataFrame:
        Select which data file to work with and load the atoms data.

    _get_atoms() -> dict[str, pd.DataFrame]:
        Get all the atoms for each residue.

    _get_residues_atoms(residues: list[str]) -> dict[str, pd.DataFrame]:
        Return a dictionary of all the residues with their atoms
        information.

    _write_msg(log: logger.logging.Logger) -> None:
        Write and log messages.

Note:
    - This class is intended to be used with gro files.
    - The script contains various methods to analyze and manipulate
        data related to residues and atoms in the system.
    - The 'param' attribute is populated with parameters from an input
        file to control various aspects of the analysis.

Example:
    data = \
        ProcessData("example.pdb", log=logger.setup_logger("update.log"))
"""


import sys
import typing
import numpy as np
import pandas as pd
import logger
import cpuconfig
import gro_to_df as grof
import read_param as param
from colors_text import TextColor as bcolors


class ProcessData:
    """
    Process the data and extract relevant sections of data for differ-
    ent residues or residue groups.

    The purpose of this script is to divide the data file and extract
    the relevant section of data. It creates separate data frames for
    different residues or groups of residues. The data is accessed th-
    rough pdb_todf.py.
    """

    info_msg: str = 'Message:\n'  # Message to pass for logging and writing
    atoms: pd.DataFrame  # All atoms dataframe
    param: dict[str, typing.Any]  # All the parameters from input file
    residues_atoms: dict[str, pd.DataFrame]  # Atoms info for each residue
    np_diameter: np.float64  # Diameter of NP, based on APTES positions
    np_depth: np.float64  # Lowest point of the np
    title: str  # Name of the system; if the file is gro
    pbc_box: str  # PBC of the system; if the file is gro

    def __init__(self,
                 fname: str,  # Name of the pdb file
                 log: logger.logging.Logger
                 ) -> None:
        """
        Initialize the ProcessData object.

        Parameters:
            fname (str): Name of the gro file.
            log (Logger): The logger object to log messages.
        """
        # Read parameters from the param file
        self.param = param.ReadParam(log=log).param

        # Get the number of cores on the host
        self.core_nr: int = self.get_nr_cores(log)

        # Load atoms data from the specified file
        self.atoms = self.__get_data(fname, log)

        # Extract atoms data for each residue and store them in a dictionary
        self.residues_atoms = self.__get_atoms()

        # Get the diameter of the NP
        self.np_diameter = self.calculate_maximum_np_radius()
        self.np_depth = self.get_np_depth()

        # Write and log the initial message
        self.__write_msg(log)

        # Empty the message
        self.info_msg = ''

    def get_nr_cores(self,
                     log: logger.logging.Logger
                     ) -> int:
        """get the number of the available cores"""
        cpu_info = cpuconfig.ConfigCpuNr(log=log)
        return cpu_info.core_nr

    def __get_data(self,
                   fname: str,  # Name of the pdb file
                   log: logger.logging.Logger
                   ) -> typing.Any:
        """
        Select which datafile to work with and load the atoms data.

        Parameters:
            fname (str): Name of the pdb file.
            log (Logger): The logger object to log messages.

        Returns:
            typing.Any: The atoms data.
        """
        if (fin := self.param['FILE']) == 'GRO':
            # Load atoms data from GRO file
            gro = grof.ReadGro(fname, log)
            atoms = gro.gro_data
            self.title = gro.title
            self.pbc_box = gro.pbc_box
        else:
            log.error(msg := f'\nFile type is not correct: {fin}!\n')
            sys.exit(f'{bcolors.FAIL}{msg}{bcolors.ENDC}')
        return atoms

    @staticmethod
    def process_chunk(chunk: np.ndarray,  # Chunk of a APTES indices
                      df_apt: pd.DataFrame  # For the APTES at the interface
                      ) -> list[int]:
        """
        Process a chunk of APTES residues to find unprotonated chains.

        Parameters:
            chunk (np.ndarray): A chunk of APTES indices to process.
            df_apt (pd.DataFrame): DataFrame containing APTES atom data.

        Returns:
            List[int]: A list of integers representing the indices of
            unprotonated APTES residues within the chunk.
        """
        # Initialize an empty list to store unprotonated APTES indices
        # in the chunk
        unprotonated_aptes_chunk: list[int] = []

        # Iterate over the APTES indices in the chunk
        for aptes_index in chunk:
            # Filter the DataFrame for the current APTES index
            df_i = df_apt[df_apt['residue_number'] == aptes_index]
            # Check if 'HN3' is present in 'atom_name' for the current
            # APTES residue
            if df_i[df_i['atom_name'].isin(['HN3'])].empty:
                # If 'HN3' is not present, add the index to the list
                # of unprotonated APTES
                unprotonated_aptes_chunk.append(aptes_index)

        # Return the list of unprotonated APTES indices in the chunk
        return unprotonated_aptes_chunk

    def __get_atoms(self) -> dict[str, pd.DataFrame]:
        """Get all the atoms for each residue.

        Returns:
            Dict[str, pd.DataFrame]: A dictionary containing pandas
            DataFrames for each residue's atoms.
        """
        # Get the names of all residues
        residues: list[str] = self.get_unique_residue_names()

        # Get the atoms for each residue and store them in a dictionary
        residues_atoms: dict[str, pd.DataFrame] = \
            self.__get_residues_atoms(residues)

        return residues_atoms

    def __get_residues_atoms(self,
                             residues: list[str]  # Name of the residues
                             ) -> dict[str, pd.DataFrame]:
        """
        Return a dictionary of all the residues with their atoms
        information.

        Parameters:
            residues (list[str]): Names of the residues.

        Returns:
            Dict[str, pd.DataFrame]: A dictionary containing pandas
            DataFrames for each residue's atoms.
        """
        residues_atoms: dict[str, pd.DataFrame] = {}  # All the atoms data
        for res in residues:
            # Filter the atoms DataFrame to get atoms belonging to each
            # residue and store them in the dictionary.
            residues_atoms[res] = self.atoms[self.atoms['residue_name'] == res]
        return residues_atoms

    def calculate_maximum_np_radius(self) -> np.float64:
        """get the maximum radius of NP, since APTES are most outward,
        here only looking at APTES residues"""
        np_diameters: list[np.float64] = []
        for aptes in self.param['aptes']:
            aptes_atoms: pd.DataFrame = self.residues_atoms[aptes]
            diameter: list[float] = []  # Save the diameters in each direction
            for xyz in ['x', 'y', 'z']:
                diameter.append(
                    aptes_atoms[xyz].max() - aptes_atoms[xyz].min())
            np_diameters.append(np.max(diameter))
        max_diameter: np.float64 = np.max(np_diameters)
        self.info_msg += \
            f'\tMaximum radius of between all NPs: `{max_diameter/2}`\n'
        return max_diameter

    def get_np_depth(self) -> np.float64:
        """returns the lowest point of the np in the water phase"""
        return np.min(self.residues_atoms['APT']['z'])

    def get_unique_residue_names(self) -> list[str]:
        """
        Get the list of the residues in the system.

        Returns:
            List[str]: A list containing the names of the residues in
            the system.
        """
        # Get the unique residue names from the 'residue_name' column
        # in the atoms DataFrame.
        residues: list[str] = list(set(self.atoms['residue_name']))
        return residues

    def __write_msg(self,
                    log: logger.logging.Logger
                    ) -> None:
        """
        Write and log messages.

        Parameters:
            log (Logger): The logger object to log the messages.
        """
        print(f'{bcolors.OKCYAN}{ProcessData.__module__}:\n'
              f'\t{self.info_msg}{bcolors.ENDC}')
        log.info(self.info_msg)


if __name__ == '__main__':
    data = ProcessData(sys.argv[1], log=logger.setup_logger('get_data.log'))

"""
write the input for plumed to se walls for ions and nanoparticles
"""

import sys
import typing
import pandas as pd

import logger
from colors_text import TextColor as bcolors


class WritePlumedInput:
    """write the input"""
    def __init__(self,
                 np_df: pd.DataFrame,  # Dataframe of nanoparticle
                 ions_df: pd.DataFrame,  # Dataframe of ions
                 log: logger.logging.Logger,
                 fout: str = 'plumed_wall.dat'  # Output file name
                 ) -> None:
        self.np_df: pd.DataFrame = np_df
        self.ions_df: pd.DataFrame = ions_df
        self.fout: str = fout
        self.write_inputs(log)

    def write_inputs(self,
                     log: logger.logging.Logger
                     ) -> None:
        """write down based on the format"""
        with open(self.fout, 'r', encoding='utf8') as f_w:
            f_w.write("# Set walls for the nanoparticles and ions\n\n")
            self.write_np_wall(f_w, log)

    def write_np_wall(self,
                      f_w: typing.IO,
                      log: logger.logging.Logger
                      ) -> None:
        """write the wall for the nanoparticles"""
        pass


if __name__ == "__main__":
    pass

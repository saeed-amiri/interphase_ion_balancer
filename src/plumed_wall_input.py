"""
write the input for plumed to se walls for ions and nanoparticles
"""

import typing

import logger
import my_tools
from colors_text import TextColor as bcolors


class WritePlumedInput:
    """write the input"""
    def __init__(self,
                 log: logger.logging.Logger,
                 fout: str = 'plumed_wall.dat'  # Output file name
                 ) -> None:
        self.fout: str = fout
        self.write_inputs(log)

    def write_inputs(self,
                     log: logger.logging.Logger
                     ) -> None:
        """write down based on the format"""
        with open(self.fout, 'w', encoding='utf8') as f_w:
            f_w.write("# Set walls for the nanoparticles and ions\n\n")
            self.write_np_wall(f_w, log)

    def write_np_wall(self,
                      f_w: typing.IO,
                      log: logger.logging.Logger
                      ) -> None:
        """write the wall for the nanoparticles"""
        self.get_atom_indices_from_index_file('index.ndx', 'CLA')

    def get_atom_indices_from_index_file(self,
                                         index_file_path: str,
                                         residue_name: str
                                         ) -> list[int]:
        """get the index of the ions from index file"""
        atom_indices: list[int] = []
        in_residue_section: bool = False
        with open(index_file_path, 'r', encoding='utf8') as index_file:
            for line in index_file:
                line = line.strip()
                if line.startswith('[') and line.endswith(']'):
                    section_name: typing.Union[str, None] = \
                        my_tools.extract_text_between_square_brackets(line)
                    if section_name is not None:
                        section_name = section_name.strip()
                        in_residue_section = section_name == residue_name
                elif in_residue_section:
                    atom_indices.extend(map(int, line.split()))

        return atom_indices


if __name__ == "__main__":
    WritePlumedInput(log=logger.setup_logger('plumed.log'))

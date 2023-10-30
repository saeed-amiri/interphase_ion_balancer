"""
write the input for restraints to ions
"""

import typing

import logger
import my_tools
from colors_text import TextColor as bcolors


class WriteRestraints:
    """write restraints fro gmx"""

    info_msg: str = 'Messages from WriteRestraints:\n'

    def __init__(self,
                 log: logger.logging.Logger,
                 parameters: dict[str, typing.Any]
                 ) -> None:
        self.write_restraints(parameters)
        self.info_msg += f'\tUsing index file: `{parameters["INDEX"]}`\n'
        self.info_msg += \
            f'\tThe restraints file is `{parameters["RESFILE"]}`\n'
        self.write_log_msg(log)

    def write_restraints(self,
                         parameters: dict[str, typing.Any]
                         ) -> None:
        """write the include file for the gromacs """
        function_id: int = 1
        ions_ndx: list[int] = \
            self.get_atom_indices_from_index_file(parameters['INDEX'],
                                                  parameters['ION'])
        with open(parameters['RESFILE'], 'w', encoding='utf8') as f_w:
            f_w.write('; Restraints on the ions\n\n')
            f_w.write('[position_restraints] \n\n')
            f_w.write('; atomnr funct fx fy fz\n')
            for atom in ions_ndx:
                f_w.write(f'{atom} {function_id} {parameters["FX"]} '
                          f'{parameters["FY"]} {parameters["FZ"]} ; '
                          f'{parameters["ION"]}\n')

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
        self.info_msg += \
            f'\tThe number ions: `{residue_name}` is `{len(atom_indices)}`\n'
        return atom_indices

    def write_log_msg(self,
                      log: logger.logging.Logger  # Name of the output file
                      ) -> None:
        """writing and logging messages from methods"""
        log.info(self.info_msg)
        print(f'{bcolors.OKBLUE}{WriteRestraints.__module__}:\n'
              f'\t{self.info_msg}\n{bcolors.ENDC}')


if __name__ == "__main__":
    import read_param as param
    LOG = logger.setup_logger('restraint.log')
    PARAM = param.ReadParam(log=LOG)
    WriteRestraints(log=LOG, parameters=PARAM.param)

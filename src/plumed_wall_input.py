"""
write the input for plumed to se walls for ions and nanoparticles
"""

import typing

import logger
import my_tools
from colors_text import TextColor as bcolors


class WritePlumedInput:
    """write the input"""

    info_msg: str = 'Messages from WritePlumedInput:\n'

    def __init__(self,
                 log: logger.logging.Logger,
                 parameters: dict[str, typing.Any],
                 ) -> None:
        self.write_inputs(log, parameters)
        self.write_log_msg(log)

    def write_inputs(self,
                     log: logger.logging.Logger,
                     parameters: dict[str, typing.Any]
                     ) -> None:
        """write down based on the format"""
        ions_ndx: list[int] = \
            self.get_atom_indices_from_index_file(parameters['INDEX'],
                                                  parameters['ION'])
        with open(parameters['OUTNAME'], 'w', encoding='utf8') as f_w:
            f_w.write("# Set walls for the nanoparticles and ions\n\n")
            f_w.write('# Walls on the nanopartcile\n')
            self.write_np_wall(f_w, parameters)
            f_w.write('# Walls on the ions\n')
            self.write_ions_wall(f_w, ions_ndx, parameters)

    def write_np_wall(self,
                      f_w: typing.IO,
                      parameters: dict[str, typing.Any]
                      ) -> None:
        """write the wall for the nanoparticles"""
        self.info_msg += f'\tThe plumed file is: `{parameters["OUTNAME"]}`\n'
        self.info_msg += f'\tThe index file is: `{parameters["INDEX"]}`\n'
        f_w.write(f'nanop: GROUP NDX_FILE={parameters["INDEX"]} '
                  f'NDX_GROUP={parameters["NP"]}\n')
        f_w.write(f'cnp: CENTER ATOMS=nanop\n'
                  f'posnp: POSITION ATOM=cnp\n'
                  f'dabs: CUSTOM ARG=posnp.z FUNC=sqrt(x*x) PERIODIC=NO\n')
        f_w.write(f'uwall_cnp: UPPER_WALLS ARG=dabs '
                  f'AT={parameters["NPCOMZ"]+0.1:.3f} '
                  f'KAPPA={parameters["KAPPA"]} '
                  f'EXP={parameters["EXP"]} EPS={parameters["EPS"]} '
                  f'OFFSET={parameters["OFFSET"]}\n')
        f_w.write(f'dwall_cnp: LOWER_WALLS ARG=dabs '
                  f'AT={parameters["NPCOMZ"]-0.05:.3f} '
                  f'KAPPA={parameters["KAPPA"]} '
                  f'EXP={parameters["EXP"]} EPS={parameters["EPS"]} '
                  f'OFFSET={parameters["OFFSET"]}\n\n')

    def write_ions_wall(self,
                        f_w: typing.IO,
                        ions_ndx: list[int],
                        parameters: dict[str, typing.Any]
                        ) -> None:
        """write the walls for the ions"""
        for i, index in enumerate(ions_ndx, start=1):
            ion_name: str = f'ions_{i}'
            ion_z_bound: tuple[float, float] = (parameters["IONLOW"] - 0.1,
                                                parameters["IONHIGH"] + 0.1)
            walls: tuple[str, str] = (f'lwall_{i}: LOWER_WALLS',
                                      f'uwall_{i}: UPPER_WALLS')
            f_w.write(f'{ion_name}: '
                      f'POSITION ATOM={index}\n')
            for wall, bound in zip(walls, ion_z_bound):
                f_w.write(f'{wall} '
                          f'ARG={ion_name}.z '
                          f'AT={bound:.3f} '
                          f'KAPPA={parameters["KAPPA"]} '
                          f'EXP={parameters["EXP"]} '
                          f'EPS={parameters["EPS"]} '
                          f'OFFSET={parameters["OFFSET"]}\n')
            f_w.write('\n')

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
        print(f'{bcolors.OKBLUE}{WritePlumedInput.__module__}:\n'
              f'\t{self.info_msg}\n{bcolors.ENDC}')


if __name__ == "__main__":
    import read_param as param
    LOG = logger.setup_logger('plumed.log')
    PARAM = param.ReadParam(log=LOG)
    WritePlumedInput(log=LOG, parameters=PARAM.param)

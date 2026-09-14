from dataclasses import dataclass
from typing import Dict, List

@dataclass
class EmissionLine:
    name: str
    rest_wave: float  # Rest frame wavelength in Angstroms
    is_main: bool = False

class LinesDB:
    """
    Database of standard emission lines and linesets.
    Replaces the kubeviz_linesdb IDL common block.
    """
    def __init__(self):
        self.lines: Dict[str, EmissionLine] = {
            'Halpha': EmissionLine('Halpha', 6562.8, is_main=True),
            'NII_6548': EmissionLine('NII_6548', 6548.05),
            'NII_6583': EmissionLine('NII_6583', 6583.45),
            'OIII_5007': EmissionLine('OIII_5007', 5006.84, is_main=True),
            'OIII_4959': EmissionLine('OIII_4959', 4958.91),
            'Hbeta': EmissionLine('Hbeta', 4861.33, is_main=True),
            'OII_3727': EmissionLine('OII_3727', 3727.09, is_main=True),
            'OII_3729': EmissionLine('OII_3729', 3729.88),
            'SII_6716': EmissionLine('SII_6716', 6716.44, is_main=True),
            'SII_6731': EmissionLine('SII_6731', 6730.82),
            'Lya': EmissionLine('Lya', 1215.67, is_main=True),
        }
        
        self.linesets: Dict[int, List[str]] = {
            1: ['Halpha', 'NII_6548', 'NII_6583'],
            2: ['OIII_5007', 'OIII_4959'],
            3: ['OII_3727', 'OII_3729'],
            4: ['Hbeta'],
            5: ['SII_6716', 'SII_6731'],
            9: ['Lya']
        }
        
    def get_lineset(self, lineset_id: int) -> List[EmissionLine]:
        if lineset_id not in self.linesets:
            return []
        return [self.lines[name] for name in self.linesets[lineset_id]]

# Singleton instance
LINES_DB = LinesDB()

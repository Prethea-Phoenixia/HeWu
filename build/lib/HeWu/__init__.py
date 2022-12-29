"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

Quick and dirty convenient functions, using the best estimate from various models
implemented in this program.
"""

from HeWu.modelBLAST1984 import airburst as _BLAST_airburst
from HeWu.modelAWG1980 import _I as _AWG_I
from HeWu.modelAWG1980 import _D_up as _AWG_D_up

__all__ = ["modelBLAST1984", "modelAWG1980", "modelBrode1970", "modelWE1984"]


def airburst(Y, HOB, GR):
    (
        PAIR,
        QAIR,
        TAAIR,
        IPTOTAL,
        DPP,
        limit1,
        XM,
        HTP,
        limit2,
        IQTOTAL,
        DPQ,
        limit3,
    ) = _BLAST_airburst(Y, HOB, GR)

    if GR < XM:
        """
        in this case, we are not in the Mach reflection region. Fall back to the
        earlier Brode model for dynamic pressure horizontal impulse component and
        positive phase duration.
        """
        IQTOTAL = _AWG_I(groundRange / 304.8, height / 304.8, Y) * 6894.76
        DPQ = _AWG_D_up(groundRange / 304.8, height / 304.8, Y) / 1e3

    return (
        PAIR,
        QAIR,
        TAAIR,
        IPTOTAL,
        DPP,
        limit1,
        XM,
        HTP,
        limit2,
        IQTOTAL,
        DPQ,
        limit3,
    )

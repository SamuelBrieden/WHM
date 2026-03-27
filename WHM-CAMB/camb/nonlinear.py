from .baseconfig import F2003Class, fortran_class
from ctypes import c_int, c_double


class NonLinearModel(F2003Class):
    """
    Abstract base class for non-linear correction models
    """
    _fields_ = [("Min_kh_nonlinear", c_double, "minimum k/h at which to apply non-linear corrections")]


halofit_original = 'original'
halofit_bird = 'bird'
halofit_peacock = 'peacock'
halofit_takahashi = 'takahashi'
halofit_mead = 'mead'
halofit_halomodel = 'halomodel'
halofit_casarini = 'casarini'
halofit_mead2015 = 'mead2015'
halofit_mead2016 = 'mead2016'
halofit_mead2020 = 'mead2020'
halofit_mead2020_feedback = 'mead2020_feedback'
halofit_brieden2023 = 'brieden2023'
halofit_brieden2023_feedback = 'brieden2023_feedback'
halofit_brieden2023_cross = 'brieden2023_cross'
halofit_brieden2023_tweaked = 'brieden2023_tweaked'
halofit_brieden2023_cross_tweaked = 'brieden2023_cross_tweaked'
halofit_halomodel_tweaked = 'halomodel_tweaked'
halofit_brieden2023_halo = 'brieden2023_halo'
halofit_brieden2023_fila = 'brieden2023_fila'
halofit_brieden2023_sheet = 'brieden2023_sheet'
halofit_brieden2023_halosphere = 'brieden2023_halosphere'
halofit_brieden2023_filasphere = 'brieden2023_filasphere'
halofit_brieden2023_sheetsphere = 'brieden2023_sheetsphere'
halofit_brieden2023_halominussphere = 'brieden2023_halominussphere'
halofit_brieden2023_filaminussphere = 'brieden2023_filaminussphere'
halofit_brieden2023_sheetminussphere = 'brieden2023_sheetminussphere'
halofit_brieden2023_halocross = 'brieden2023_halocross'
halofit_brieden2023_filacross = 'brieden2023_filacross'
halofit_brieden2023_sheetcross = 'brieden2023_sheetcross'
halofit_brieden2023_halospherecross = 'brieden2023_halospherecross'
halofit_brieden2023_filaspherecross = 'brieden2023_filaspherecross'
halofit_brieden2023_sheetspherecross = 'brieden2023_sheetspherecross'
halofit_brieden2023_cross_profilecusp = 'brieden2023_cross_profilecored'
halofit_brieden2023_cross_profileconc = 'brieden2023_cross_profilecusp'
halofit_brieden2023_cross_nobias = 'brieden2023_cross_nobias'
halofit_brieden2023_cross_nohalobias = 'brieden2023_cross_nohalobias'
halofit_brieden2023_cross_tinkerhalobias = 'brieden2023_cross_tinkerhalobias'
halofit_brieden2023_cross_fnumasscutlow = 'brieden2023_cross_fnumasscutlow'
halofit_brieden2023_cross_fnumasscuthigh = 'brieden2023_cross_fnumasscuthigh'

halofit_default = halofit_mead2020

halofit_version_names = {halofit_original: 1,
                         halofit_bird: 2,
                         halofit_peacock: 3,
                         halofit_takahashi: 4,
                         halofit_mead: 5,
                         halofit_halomodel: 6,
                         halofit_casarini: 7,
                         halofit_mead2015: 8,
                         halofit_mead2016: 5,
                         halofit_mead2020: 9,
                         halofit_mead2020_feedback: 10,
                         halofit_brieden2023: 11,
                         halofit_brieden2023_feedback: 12,
                         halofit_brieden2023_cross: 13,
                         halofit_brieden2023_tweaked: 14,
                         halofit_brieden2023_cross_tweaked: 15,
                         halofit_halomodel_tweaked: 16,
                         halofit_brieden2023_halo: 17,
                         halofit_brieden2023_fila: 18,
                         halofit_brieden2023_sheet: 19,
                         halofit_brieden2023_halosphere: 20,
                         halofit_brieden2023_filasphere: 21,
                         halofit_brieden2023_sheetsphere: 22,
                         halofit_brieden2023_halominussphere: 23,
                         halofit_brieden2023_filaminussphere: 24,
                         halofit_brieden2023_sheetminussphere: 25,
                         halofit_brieden2023_halocross: 26,
                         halofit_brieden2023_filacross: 27,
                         halofit_brieden2023_sheetcross: 28,
                         halofit_brieden2023_halospherecross: 29,
                         halofit_brieden2023_filaspherecross: 30,
                         halofit_brieden2023_sheetspherecross: 31,
                         halofit_brieden2023_cross_profilecusp: 32,
                         halofit_brieden2023_cross_profileconc: 33,
                         halofit_brieden2023_cross_nobias: 34,
                         halofit_brieden2023_cross_nohalobias: 35,
                         halofit_brieden2023_cross_tinkerhalobias: 36,
                         halofit_brieden2023_cross_fnumasscutlow: 37,
                         halofit_brieden2023_cross_fnumasscuthigh: 38,
                         
                        }


@fortran_class
class Halofit(NonLinearModel):
    """
    Various specific approximate non-linear correction models based on HaloFit.
    """
    _fields_ = [
        ("halofit_version", c_int, {"names": halofit_version_names}),
        ("HMCode_A_baryon", c_double, "HMcode parameter A_baryon"),
        ("HMCode_eta_baryon", c_double, "HMcode parameter eta_baryon"),
        ("HMCode_logT_AGN", c_double, "HMcode parameter log10(T_AGN/K)")
    ]

    _fortran_class_module_ = 'NonLinear'
    _fortran_class_name_ = 'THalofit'

    def get_halofit_version(self):
        return self.halofit_version

    def set_params(self, halofit_version=halofit_default, HMCode_A_baryon=3.13, HMCode_eta_baryon=0.603,
                   HMCode_logT_AGN=7.8):
        """
        Set the halofit model for non-linear corrections.

        :param halofit_version: One of

            - original: `astro-ph/0207664 <https://arxiv.org/abs/astro-ph/0207664>`_
            - bird: `arXiv:1109.4416 <https://arxiv.org/abs/1109.4416>`_
            - peacock: `Peacock fit <http://www.roe.ac.uk/~jap/haloes/>`_
            - takahashi: `arXiv:1208.2701 <https://arxiv.org/abs/1208.2701>`_
            - mead: HMCode `arXiv:1602.02154 <https://arxiv.org/abs/1602.02154>`_
            - halomodel: basic halomodel
            - casarini: PKequal `arXiv:0810.0190 <https://arxiv.org/abs/0810.0190>`_, `arXiv:1601.07230 <https://arxiv.org/abs/1601.07230>`_
            - mead2015: original 2015 version of HMCode `arXiv:1505.07833 <https://arxiv.org/abs/1505.07833>`_
            - mead2016: Alias for 'mead'.
            - mead2020: 2020 version of HMcode `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
            - mead2020_feedback: 2020 version of HMcode with baryonic feedback `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
        :param HMCode_A_baryon: HMcode parameter A_baryon. Default 3.13. Used only in models mead2015 and mead2016 (and its alias mead).
        :param HMCode_eta_baryon: HMcode parameter eta_baryon. Default 0.603. Used only in mead2015 and mead2016 (and its alias mead).
        :param HMCode_logT_AGN: HMcode parameter logT_AGN. Default 7.8. Used only in model mead2020_feedback.
        """
        self.halofit_version = halofit_version
        self.HMCode_A_baryon = HMCode_A_baryon
        self.HMCode_eta_baryon = HMCode_eta_baryon
        self.HMCode_logT_AGN = HMCode_logT_AGN


@fortran_class
class SecondOrderPK(NonLinearModel):
    """
    Third-order Newtonian perturbation theory results for the non-linear correction.
    Only intended for use at very high redshift (z>10) where corrections are perturbative, it will not give
    sensible results at low redshift.

    See Appendix F of `astro-ph/0702600 <https://arxiv.org/abs/astro-ph/0702600>`_ for equations and references.

    Not intended for production use, it's mainly to serve as an example alternative non-linear model implementation.
    """

    _fortran_class_module_ = 'SecondOrderPK'
    _fortran_class_name_ = 'TSecondOrderPK'

    def set_params(self):
        pass

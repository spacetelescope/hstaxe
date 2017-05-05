from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
def get_HRC_PR200L_WCS():
    """
    Defines parameters for the ACS/HRC/PR200L slitless mode

    Returns
    -------
    params: dict
        slitless mode parameters
    """
    wcs_keys = {}
    # grism image "j97b06sdq":
    # / World Coordinate System and Related Parameters
    wcs_keys['grism'] = [['WCSAXES', 2, 'number of World Coordinate System axes'],
                        ['CRPIX1',   6.520000000000E+02, 'x-coordinate of reference pixel'],
                        ['CRPIX2',   5.120000000000E+02, 'y-coordinate of reference pixel'],
                        ['CRVAL1',   3.230683853815E+02, 'first axis value at reference pixel'],
                        ['CRVAL2',   2.539875952324E-01, 'second axis value at reference pixel'],
                        ['CTYPE1', 'RA---TAN', 'the coordinate type for the first axis'],
                        ['CTYPE2', 'DEC--TAN', 'the coordinate type for the second axis'],
                        ['CD1_1',           -9.554E-07, 'partial of first axis coordinate w.r.t. x'],
                        ['CD1_2',          6.73946E-06, 'partial of first axis coordinate w.r.t. y'],
                        ['CD2_1',          7.85265E-06, 'partial of second axis coordinate w.r.t. x'],
                        ['CD2_2',          1.49907E-06, 'partial of second axis coordinate w.r.t. y'],
                        ['LTV1',        0.0000000E+00, 'offset in X to subsection start'],
                        ['LTV2',        0.0000000E+00, 'offset in Y to subsection start'],
                        ['LTM1_1',                  1.0, 'reciprocal of sampling rate in X'],
                        ['LTM2_2',                  1.0, 'reciprocal of sampling rate in Y'],
                        ['ORIENTAT',              77.4597, 'position angle of image y axis (deg. e of n)'],
                        ['RA_APER',   3.230683853815E+02, 'RA of aperture reference position'],
                        ['DEC_APER',   2.539875952324E-01, 'Declination of aperture reference position'],
                        ['PA_APER',              77.4597, 'Position Angle of reference aperture center (de'],
                        ['VAFACTOR',   1.000109352103E+00, 'velocity aberration plate scale factor'],
                        ['EXPNAME', 'j97b06sdq', 'exposure identifier']]

    wcs_keys['drizzle'] = [
            ['DRZCNUM', 15, 'Number of coefficients per coordinate'],
            ['DRZSCALE', 0.025, 'Scale for drizzling'],
            ['DRZ2X01', 2.08138254, 'Drizzle coefficient 01 in X'],
            ['DRZ2X02', 1.1315262, 'Drizzle coefficient 02 in X'],
            ['DRZ2X03', -0.0015648471, 'Drizzle coefficient 03 in X'],
            ['DRZ2X04', -3.7322242e-06, 'Drizzle coefficient 04 in X'],
            ['DRZ2X05', 1.0109418e-05, 'Drizzle coefficient 05 in X'],
            ['DRZ2X06', -8.1409948e-07, 'Drizzle coefficient 06 in X'],
            ['DRZ2X07', 1.3075234e-10, 'Drizzle coefficient 07 in X'],
            ['DRZ2X08', 2.3539512e-11, 'Drizzle coefficient 08 in X'],
            ['DRZ2X09', 4.3735033e-10, 'Drizzle coefficient 09 in X'],
            ['DRZ2X10', -1.0632533e-11, 'Drizzle coefficient 10 in X'],
            ['DRZ2X11', 2.7818193e-13, 'Drizzle coefficient 10 in X'],
            ['DRZ2X12', -2.0905397e-12, 'Drizzle coefficient 10 in X'],
            ['DRZ2X13', -1.6363199e-12, 'Drizzle coefficient 10 in X'],
            ['DRZ2X14', 1.4459954e-13, 'Drizzle coefficient 10 in X'],
            ['DRZ2X15', 4.6694459e-14, 'Drizzle coefficient 10 in X'],
            ['DRZ2Y01', -2.81487211, 'Drizzle coefficient 01 in Y'],
            ['DRZ2Y02', 0.11637744, 'Drizzle coefficient 02 in Y'],
            ['DRZ2Y03', 0.99369911, 'Drizzle coefficient 03 in Y'],
            ['DRZ2Y04', 1.5666096e-06, 'Drizzle coefficient 04 in Y'],
            ['DRZ2Y05', -1.5948952e-06, 'Drizzle coefficient 05 in Y'],
            ['DRZ2Y06', 1.1298673e-05, 'Drizzle coefficient 06 in Y'],
            ['DRZ2Y07',  4.5586465e-10, 'Drizzle coefficient 07 in Y'],
            ['DRZ2Y08',  3.300937e-10, 'Drizzle coefficient 08 in Y'],
            ['DRZ2Y09', -6.4272644e-10, 'Drizzle coefficient 09 in Y'],
            ['DRZ2Y10', 4.9637125e-10, 'Drizzle coefficient 10 in Y'],
            ['DRZ2Y11', 4.1921892e-13, 'Drizzle coefficient 10 in Y'],
            ['DRZ2Y12', -1.6327524e-12, 'Drizzle coefficient 10 in Y'],
            ['DRZ2Y13', 8.2560585e-13, 'Drizzle coefficient 10 in Y'],
            ['DRZ2Y14', -5.5250041e-13, 'Drizzle coefficient 10 in Y'],
            ['DRZ2Y15', -1.0909078e-12, 'Drizzle coefficient 10 in Y'],
            ]

    wcs_keys['direct'] = [
        # direct image "j97b06scq":
        # / World Coordinate System and Related Parameters
        ['WCSAXES', 2, 'number of World Coordinate System axes'],
        ['CRPIX1', 5.120000000000E+02, 'x-coordinate of reference pixel'],
        ['CRPIX2', 5.120000000000E+02, 'y-coordinate of reference pixel'],
        ['CRVAL1', 3.230683853809E+02, 'first axis value at reference pixel'],
        ['CRVAL2',   2.539875951692E-01, 'second axis value at reference pixel'],
        ['CTYPE1', 'RA---TAN', 'the coordinate type for the first axis'],
        ['CTYPE2', 'DEC--TAN', 'the coordinate type for the second axis'],
        ['CD1_1', -9.47228E-07, 'partial of first axis coordinate w.r.t. x'],
        ['CD1_2', 6.73225E-06, 'partial of first axis coordinate w.r.t. y'],
        ['CD2_1', 7.84049E-06, 'partial of second axis coordinate w.r.t. x'],
        ['CD2_2', 1.51604E-06, 'partial of second axis coordinate w.r.t. y'],
        ['LTV1', 0.0000000E+00, 'offset in X to subsection start'],
        ['LTV2', 0.0000000E+00, 'offset in Y to subsection start'],
        ['LTM1_1', 1.0, 'reciprocal of sampling rate in X'],
        ['LTM2_2', 1.0, 'reciprocal of sampling rate in Y'],
        ['ORIENTAT', 77.3092, 'position angle of image y axis (deg. e of n)'],
        ['RA_APER', 3.230683853809E+02, 'RA of aperture reference position'],
        ['DEC_APER', 2.539875951692E-01, 'Declination of aperture reference position'],
        ['PA_APER', 77.3092, 'Position Angle of reference aperture center (de'],
        ['VAFACTOR', 1.000110148008E+00, 'velocity aberration plate scale factor'],
        ['EXPNAME ', 'j97b06scq', 'exposure identifier'],
        ]

    wcs_keys['dimension'] = [1024, 1024]
    return wcs_keys


def get_SBC_PR110L_WCS():
    """Defines parameters for the ACS/SBC/PR110L slitless mode

    Returns
    -------
    params: dict
        slitless mode parameters
    """
    wcs_keys = {}
    # grism image "j97b11slq":
    # / World Coordinate System and Related Parameters
    wcs_keys['grism'] = [
        ['WCSAXES',                    2 ,'number of World Coordinate System axes'],
        ['CRPIX1',                512.0 ,'x-coordinate of reference pixel'],
        ['CRPIX2',                512.0 ,'y-coordinate of reference pixel'],
        ['CRVAL1',    323.0683271777468 ,'first axis value at reference pixel'],
        ['CRVAL2',   0.2552839595533656 ,'second axis value at reference pixel'],
        ['CTYPE1', 'RA---TAN'           ,'the coordinate type for the first axis'],
        ['CTYPE2', 'DEC--TAN'           ,'the coordinate type for the second axis'],
        ['CD1_1', -7.161023752814579E-06 ,'partial of first axis coordinate w.r.t. x'],
        ['CD1_2', 4.749202573663841E-06 ,'partial of first axis coordinate w.r.t. y'],
        ['CD2_1', 6.047100892329302E-06 ,'partial of second axis coordinate w.r.t. x'],
        ['CD2_2', 6.879329383769482E-06 ,'partial of second axis coordinate w.r.t. y'],
        ['LTV1',                  0.0 ,'offset in X to subsection start'],
        ['LTV2',                  0.0 ,'offset in Y to subsection start'],
        ['LTM1_1',                  1.0 ,'reciprocal of sampling rate in X'],
        ['LTM2_2',                  1.0 ,'reciprocal of sampling rate in Y'],
        ['ORIENTAT',    34.63527041715778 ,'position angle of image y axis (deg. e of n)'],
        ['RA_APER',   3.230684171774E+02 ,'RA of aperture reference position'],
        ['DEC_APER',   2.539909262239E-01 ,'Declination of aperture reference position'],
        ['PA_APER',              34.6353 ,'Position Angle of reference aperture center (de'],
        ['VAFACTOR',   1.000053103518E+00 ,'velocity aberration plate scale factor'],
        ['EXPNAME', 'j97b11slq' ,'exposure identifier'],
        ]

    wcs_keys['drizzle'] = [
        ['DRZCNUM',15,'Number of coefficients per coordinate'],
        ['DRZSCALE',0.025,'Scale for drizzling'],
        ['DRZ2X01', 1.04924098,'Drizzle coefficient 01 in X'],
        ['DRZ2X02', 1.3453513,'Drizzle coefficient 02 in X'],
        ['DRZ2X03', 0.020481231,'Drizzle coefficient 03 in X'],
        ['DRZ2X04', -1.844869e-05,'Drizzle coefficient 04 in X'],
        ['DRZ2X05', 1.6256479e-05,'Drizzle coefficient 05 in X'],
        ['DRZ2X06', 5.8220551e-06,'Drizzle coefficient 06 in X'],
        ['DRZ2X07', 4.0145121e-09,'Drizzle coefficient 07 in X'],
        ['DRZ2X08', 2.3509428e-09,'Drizzle coefficient 08 in X'],
        ['DRZ2X09', 2.3509428e-09,'Drizzle coefficient 09 in X'],
        ['DRZ2X10', -7.6725898e-09,'Drizzle coefficient 10 in X'],
        ['DRZ2X11', 1.4285034e-11,'Drizzle coefficient 10 in X'],
        ['DRZ2X12', 1.6199441e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X13', 2.1826236e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X14', -7.9748938e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X15', -7.9748938e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2Y01', -4.29830391,'Drizzle coefficient 01 in Y'],
        ['DRZ2Y02', 0.10788573,'Drizzle coefficient 02 in Y'],
        ['DRZ2Y03', 1.2035841,'Drizzle coefficient 03 in Y'],
        ['DRZ2Y04', -1.4603854e-06,'Drizzle coefficient 04 in Y'],
        ['DRZ2Y05', 2.437245e-06,'Drizzle coefficient 05 in Y'],
        ['DRZ2Y06', 1.2681371e-05,'Drizzle coefficient 06 in Y'],
        ['DRZ2Y07', 1.5990346e-09,'Drizzle coefficient 07 in Y'],
        ['DRZ2Y08', -8.8217959e-09,'Drizzle coefficient 08 in Y'],
        ['DRZ2Y09', -1.0642736e-09,'Drizzle coefficient 09 in Y'],
        ['DRZ2Y10', 1.6438327e-09,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y11', 3.5539949e-12,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y12', -9.4676107e-13,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y13', 1.8894413e-12,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y14', 8.9042909e-12,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y15', -2.9544276e-12,'Drizzle coefficient 10 in Y'],
        ]

    #direct image "j97b11sjq":
    #/ World Coordinate System and Related Parameters
    wcs_keys['direct'] = [
        ['WCSAXES',                    2 ,'number of World Coordinate System axes'],
        ['CRPIX1',                512.0 ,'x-coordinate of reference pixel'],
        ['CRPIX2',                512.0 ,'y-coordinate of reference pixel'],
        ['CRVAL1',    323.0689488179553 ,'first axis value at reference pixel'],
        ['CRVAL2',    0.254758917978487 ,'second axis value at reference pixel'],
        ['CTYPE1', 'RA---TAN'           ,'the coordinate type for the first axis'],
        ['CTYPE2', 'DEC--TAN'           ,'the coordinate type for the second axis'],
        ['CD1_1', -7.161034143224408E-06 ,'partial of first axis coordinate w.r.t. x'],
        ['CD1_2', 4.749209388017996E-06 ,'partial of first axis coordinate w.r.t. y'],
        ['CD2_1', 6.04710961189569E-06 ,'partial of second axis coordinate w.r.t. x'],
        ['CD2_2', 6.879339303408634E-06 ,'partial of second axis coordinate w.r.t. y'],
        ['LTV1',                  0.0 ,'offset in X to subsection start'],
        ['LTV2',                  0.0 ,'offset in Y to subsection start'],
        ['LTM1_1',                  1.0 ,'reciprocal of sampling rate in X'],
        ['LTM2_2',                  1.0 ,'reciprocal of sampling rate in Y'],
        ['ORIENTAT',    34.69259860641057 ,'position angle of image y axis (deg. e of n)'],
        ['RA_APER',   3.230684171767E+02 ,'RA of aperture reference position'],
        ['DEC_APER',   2.539909261495E-01 ,'Declination of aperture reference position'],
        ['PA_APER',              34.6926 ,'Position Angle of reference aperture center (de'],
        ['VAFACTOR',   1.000054545527E+00 ,'velocity aberration plate scale factor'],
        ['EXPNAME', 'j97b11sjq' ,'exposure identifier'],
        ]
    wcs_keys['dimension'] = [1024, 1024]
    return wcs_keys

def get_SBC_PR130L_WCS():
    """
    Defines parameters for the ACS/SBC/PR130L slitless mode

    @return: slitless mode parameters
    @rtype: dictionary
    """
    wcs_keys = {}
    #grism image "j97b11smq":
    #/ World Coordinate System and Related Parameters
    wcs_keys['grism'] = [
        ['WCSAXES ',                    2 ,'number of World Coordinate System axes'],
        ['CRPIX1',                512.0 ,'x-coordinate of reference pixel'],
        ['CRPIX2',                512.0 ,'y-coordinate of reference pixel'],
        ['CRVAL1',    323.0683271781467 ,'first axis value at reference pixel'],
        ['CRVAL2',   0.2552839595928655 ,'second axis value at reference pixel'],
        ['CTYPE1', 'RA---TAN'           ,'the coordinate type for the first axis'],
        ['CTYPE2', 'DEC--TAN'           ,'the coordinate type for the second axis'],
        ['CD1_1', -7.161016148187144E-06 ,'partial of first axis coordinate w.r.t. x'],
        ['CD1_2', 4.749197435716537E-06 ,'partial of first axis coordinate w.r.t. y'],
        ['CD2_1', 6.047094422625643E-06 ,'partial of second axis coordinate w.r.t. x'],
        ['CD2_2', 6.879322023732403E-06 ,'partial of second axis coordinate w.r.t. y'],
        ['LTV1',                  0.0 ,'offset in X to subsection start'],
        ['LTV2',                  0.0 ,'offset in Y to subsection start'],
        ['LTM1_1',                  1.0 ,'reciprocal of sampling rate in X'],
        ['LTM2_2',                  1.0 ,'reciprocal of sampling rate in Y'],
        ['ORIENTAT',    34.63527041715778 ,'position angle of image y axis (deg. e of n)'],
        ['RA_APER',   3.230684171778E+02 ,'RA of aperture reference position '],
        ['DEC_APER',   2.539909262634E-01 ,'Declination of aperture reference position'],
        ['PA_APER',              34.6353 ,'Position Angle of reference aperture center (de'],
        ['VAFACTOR',   1.000052033576E+00 ,'velocity aberration plate scale factor'],
        ['EXPNAME', 'j97b11smq' ,'exposure identifier'],
        ]


    wcs_keys['drizzle'] = [
        ['DRZCNUM',15,'Number of coefficients per coordinate'],
        ['DRZSCALE',0.025,'Scale for drizzling'],
        ['DRZ2X01', 1.04923371,'Drizzle coefficient 01 in X'],
        ['DRZ2X02', 1.3453499,'Drizzle coefficient 02 in X'],
        ['DRZ2X03', 0.020481207,'Drizzle coefficient 03 in X'],
        ['DRZ2X04', -1.8448652e-05,'Drizzle coefficient 04 in X'],
        ['DRZ2X05', 1.6256445e-05,'Drizzle coefficient 05 in X'],
        ['DRZ2X06', 5.822046e-06,'Drizzle coefficient 06 in X'],
        ['DRZ2X07', 4.0144937e-09,'Drizzle coefficient 07 in X'],
        ['DRZ2X08', 2.3509343e-09,'Drizzle coefficient 08 in X'],
        ['DRZ2X09', -5.7698199e-09,'Drizzle coefficient 09 in X'],
        ['DRZ2X10', -7.6725457e-09,'Drizzle coefficient 10 in X'],
        ['DRZ2X11', 1.4284973e-11,'Drizzle coefficient 10 in X'],
        ['DRZ2X12', 1.6199372e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X13', 2.1826142e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X14', -7.9748596e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X15', -3.9227763e-11,'Drizzle coefficient 10 in X'],
        ['DRZ2Y01', -4.29829370,'Drizzle coefficient 01 in Y'],
        ['DRZ2Y02', 0.10788562,'Drizzle coefficient 02 in Y'],
        ['DRZ2Y03', 1.2035828,'Drizzle coefficient 03 in Y'],
        ['DRZ2Y04', -1.4603816e-06,'Drizzle coefficient 04 in Y'],
        ['DRZ2Y05', 2.4372417e-06,'Drizzle coefficient 05 in Y'],
        ['DRZ2Y06', 1.2681344e-05,'Drizzle coefficient 06 in Y'],
        ['DRZ2Y07', 1.5990282e-09,'Drizzle coefficient 07 in Y'],
        ['DRZ2Y08', -8.8217677e-09,'Drizzle coefficient 08 in Y'],
        ['DRZ2Y09', -1.0642737e-09,'Drizzle coefficient 09 in Y'],
        ['DRZ2Y10', 1.643828e-09 ,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y11', 3.5539797e-12,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y12', -9.4675702e-13,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y13', 1.8894332e-12,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y14', 8.9042528e-12,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y15', -2.954415e-12,'Drizzle coefficient 10 in Y'],
        ]

    wcs_keys['direct'] = [
        #direct image "j97b11sjq":
        #/ World Coordinate System and Related Parameters
        ['WCSAXES',2 ,'number of World Coordinate System axes'],
        ['CRPIX1',512.0 ,'x-coordinate of reference pixel'],
        ['CRPIX2',512.0 ,'y-coordinate of reference pixel'],
        ['CRVAL1',323.0689488179553 ,'first axis value at reference pixel'],
        ['CRVAL2',0.254758917978487 ,'second axis value at reference pixel'],
        ['CTYPE1','RA---TAN','the coordinate type for the first axis'],
        ['CTYPE2','DEC--TAN','the coordinate type for the second axis'],
        ['CD1_1',-7.161034143224408E-06 ,'partial of first axis coordinate w.r.t. x'],
        ['CD1_2',4.749209388017996E-06 ,'partial of first axis coordinate w.r.t. y'],
        ['CD2_1',6.04710961189569E-06 ,'partial of second axis coordinate w.r.t. x'],
        ['CD2_2',6.879339303408634E-06 ,'partial of second axis coordinate w.r.t. y'],
        ['LTV1',0.0 ,'offset in X to subsection start'],
        ['LTV2',0.0 ,'offset in Y to subsection start'],
        ['LTM1_1',1.0 ,'reciprocal of sampling rate in X'],
        ['LTM2_2',1.0 ,'reciprocal of sampling rate in Y'],
        ['ORIENTAT',34.69259860641057 ,'position angle of image y axis (deg. e of n)'],
        ['RA_APER',3.230684171767E+02 ,'RA of aperture reference position'],
        ['DEC_APER',2.539909261495E-01 ,'Declination of aperture reference position'],
        ['PA_APER',34.6926 ,'Position Angle of reference aperture center (de'],
        ['VAFACTOR',1.000054545527E+00 ,'velocity aberration plate scale factor'],
        ['EXPNAME','j97b11sjq' ,'exposure identifier'],
        ]
    wcs_keys['dimension'] = [1024, 1024]
    return wcs_keys


def get_HRC_G800L_WCS():
    """
    Defines parameters for the ACS/HRC/G800L slitless mode

    @return: slitless mode parameters
    @rtype: dictionary
    """
    wcs_keys = {}
    #grism image "j8m820leq":
    #/ World Coordinate System and Related Parameters
    wcs_keys['grism'] = [
        ['WCSAXES',2 ,'number of World Coordinate System axes'],
        ['CRPIX1',5.120000000000E+02 ,'x-coordinate of reference pixel'],
        ['CRPIX2',5.120000000000E+02 ,'y-coordinate of reference pixel'],
        ['CRVAL1',5.322491146968E+01 ,'first axis value at reference pixel'],
        ['CRVAL2',-2.782350140649E+01 ,'second axis value at reference pixel'],
        ['CTYPE1','RA---TAN','the coordinate type for the first axis'],
        ['CTYPE2','DEC--TAN','the coordinate type for the second axis'],
        ['CD1_1',6.03764E-06 ,'partial of first axis coordinate w.r.t. x'],
        ['CD1_2',4.95891E-06 ,'partial of first axis coordinate w.r.t. y'],
        ['CD2_1',5.09479E-06 ,'partial of second axis coordinate w.r.t. x'],
        ['CD2_2',-4.79903E-06 ,'partial of second axis coordinate w.r.t. y'],
        ['LTV1',0.0000000E+00 ,'offset in X to subsection start'],
        ['LTV2',0.0000000E+00 ,'offset in Y to subsection start'],
        ['LTM1_1',1.0 ,'reciprocal of sampling rate in X'],
        ['LTM2_2',1.0 ,'reciprocal of sampling rate in Y'],
        ['ORIENTAT',134.061 ,'position angle of image y axis (deg. e of n)'],
        ['RA_APER',5.322491146968E+01 ,'RA of aperture reference position'],
        ['DEC_APER',-2.782350140649E+01 ,'Declination of aperture reference position'],
        ['PA_APER',134.061 ,'Position Angle of reference aperture center (de'],
        ['VAFACTOR',1.000045790665E+00 ,'velocity aberration plate scale factor'],
        ['EXPNAME', 'j8m820leq ' ,'exposure identifier']
        ]


    wcs_keys['drizzle'] = [
        ['DRZCNUM',15,'Number of coefficients per coordinate'],
        ['DRZSCALE',0.025,'Scale for drizzling'],
        ['DRZ2X01', 1.99990688,'Drizzle coefficient 01 in X'],
        ['DRZ2X02', 1.1312952,'Drizzle coefficient 02 in X'],
        ['DRZ2X03', -0.0013998074,'Drizzle coefficient 03 in X'],
        ['DRZ2X04', -3.7547287e-06,'Drizzle coefficient 04 in X'],
        ['DRZ2X05', 1.0207302e-05,'Drizzle coefficient 05 in X'],
        ['DRZ2X06', -8.1810015e-07,'Drizzle coefficient 06 in X'],
        ['DRZ2X07', 1.0776155e-10,'Drizzle coefficient 07 in X'],
        ['DRZ2X08', -3.5401801e-11,'Drizzle coefficient 08 in X'],
        ['DRZ2X09', 4.6434155e-10,'Drizzle coefficient 09 in X'],
        ['DRZ2X10', -2.226957e-11,'Drizzle coefficient 10 in X'],
        ['DRZ2X11', 3.6592701e-13,'Drizzle coefficient 10 in X'],
        ['DRZ2X12', -2.5148205e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X13', -1.5241192e-12,'Drizzle coefficient 10 in X'],
        ['DRZ2X14', -2.0761802e-13,'Drizzle coefficient 10 in X'],
        ['DRZ2X15', 4.1691635e-14,'Drizzle coefficient 10 in X'],
        ['DRZ2Y01', -2.82194592,'Drizzle coefficient 01 in Y'],
        ['DRZ2Y02', 0.11656498,'Drizzle coefficient 02 in Y'],
        ['DRZ2Y03', 0.99377446,'Drizzle coefficient 03 in Y'],
        ['DRZ2Y04', 1.6018684e-06,'Drizzle coefficient 04 in Y'],
        ['DRZ2Y05', -1.6023185e-06,'Drizzle coefficient 05 in Y'],
        ['DRZ2Y06', 1.1354663e-05,'Drizzle coefficient 06 in Y'],
        ['DRZ2Y07', 4.5674659e-10,'Drizzle coefficient 07 in Y'],
        ['DRZ2Y08', 2.9705672e-10,'Drizzle coefficient 08 in Y'],
        ['DRZ2Y09', -5.8000163e-10,'Drizzle coefficient 09 in Y'],
        ['DRZ2Y10', 7.2709968e-10,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y11', 3.0293949e-13,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y12', -1.5578453e-12,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y13', 8.3043208e-13,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y14', -5.5750209e-13,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y15', -1.2661919e-12,'Drizzle coefficient 10 in Y'],
        ]

    wcs_keys['direct'] = None
    wcs_keys['dimension'] = [1024, 1024]
    return wcs_keys

def get_WFC_G800L_WCS():
    """
    Defines parameters for the ACS/WFC/G800L slitless mode

    @return: slitless mode parameters
    @rtype: dictionary
    """
    wcs_keys = {}
    #grism image "j8qq12kgq":
    #/ World Coordinate System and Related Parameters
    wcs_keys['grism'] = [
        ['WCSAXES ',2,'number of World Coordinate System axes'],
        ['CRPIX1',2.048000000000E+03 ,'x-coordinate of reference pixel'],
        ['CRPIX2 ',1.024000000000E+03 ,'y-coordinate of reference pixel'],
        ['CRVAL1',5.317727241178E+01 ,'first axis value at reference pixel'],
        ['CRVAL2',-2.779882660148E+01 ,'second axis value at reference pixel'],
        ['CTYPE1','RA---TAN','the coordinate type for the first axis'],
        ['CTYPE2','DEC--TAN','the coordinate type for the second axis'],
        ['CD1_1',-8.4686E-06 ,'partial of first axis coordinate w.r.t. x'],
        ['CD1_2',-1.16183E-05 ,'partial of first axis coordinate w.r.t. y'],
        ['CD2_1',-1.09568E-05 ,'partial of second axis coordinate w.r.t. x'],
        ['CD2_2',7.76629E-06 ,'partial of second axis coordinate w.r.t. y'],
        ['LTV1',0.0000000E+00 ,'offset in X to subsection start'],
        ['LTV2',0.0000000E+00 ,'offset in Y to subsection start'],
        ['LTM1_1',1.0 ,'reciprocal of sampling rate in X'],
        ['LTM2_2',1.0 ,'reciprocal of sampling rate in Y'],
        ['ORIENTAT',-56.2392 ,'position angle of image y axis (deg. e of n)'],
        ['RA_APER',5.316373797542E+01 ,'RA of aperture reference position'],
        ['DEC_APER',-2.779155504028E+01 ,'Declination of aperture reference position'],
        ['PA_APER',-56.4791 ,'Position Angle of reference aperture center (de'],
        ['VAFACTOR',1.000027253644E+00 ,'velocity aberration plate scale factor'],
        ['EXPNAME','j8qq12kgq' ,'exposure identifier'],
        ]

    wcs_keys['drizzle'] = [
        ['DRZCNUM',15,'Number of coefficients per coordinate'],
        ['DRZSCALE',0.05,'Scale for drizzling'],
        ['DRZ2X01', 29.20952271 ,'Drizzle coefficient 01 in X'],
        ['DRZ2X02', 0.98463856,'Drizzle coefficient 02 in X'],
        ['DRZ2X03', 0.047121902,'Drizzle coefficient 03 in X'],
        ['DRZ2X04', 8.2405479e-06,'Drizzle coefficient 04 in X'],
        ['DRZ2X05', -7.1109122e-06,'Drizzle coefficient 05 in X'],
        ['DRZ2X06', 1.7714826e-06,'Drizzle coefficient 06 in X'],
        ['DRZ2X07', -4.6293307e-10,'Drizzle coefficient 07 in X'],
        ['DRZ2X08', -1.243901e-10,'Drizzle coefficient 08 in X'],
        ['DRZ2X09', -5.3285875e-10,'Drizzle coefficient 09 in X'],
        ['DRZ2X10', 5.1490811e-11,'Drizzle coefficient 10 in X'],
        ['DRZ2X11', 1.6734254e-14,'Drizzle coefficient 11 in X'],
        ['DRZ2X12', 3.425828e-14,'Drizzle coefficient 12 in X'],
        ['DRZ2X13', 9.5688062e-14 ,'Drizzle coefficient 13 in X'],
        ['DRZ2X14', -1.6229259e-14,'Drizzle coefficient 14 in X'],
        ['DRZ2X15', 1.2711148e-13,'Drizzle coefficient 15 in X'],
        ['DRZ2Y01', 1047.90670925,'Drizzle coefficient 01 in Y'],
        ['DRZ2Y02', 0.040785422,'Drizzle coefficient 02 in Y'],
        ['DRZ2Y03', 0.97161774,'Drizzle coefficient 03 in Y'],
        ['DRZ2Y04', -2.5332551e-06,'Drizzle coefficient 04 in Y'],
        ['DRZ2Y05', 5.9183197e-06,'Drizzle coefficient 05 in Y'],
        ['DRZ2Y06', -9.4306843e-06,'Drizzle coefficient 06 in Y'],
        ['DRZ2Y07', 7.3674246e-11,'Drizzle coefficient 07 in Y'],
        ['DRZ2Y08', -4.3916951e-10,'Drizzle coefficient 08 in Y'],
        ['DRZ2Y09', -5.371583e-11,'Drizzle coefficient 09 in Y'],
        ['DRZ2Y10', -3.8747876e-10,'Drizzle coefficient 10 in Y'],
        ['DRZ2Y11', -1.4892746e-14,'Drizzle coefficient 11 in Y'],
        ['DRZ2Y12', -3.1028203e-14,'Drizzle coefficient 12 in Y'],
        ['DRZ2Y13', -1.024679e-13,'Drizzle coefficient 13 in Y'],
        ['DRZ2Y14', 2.9690206e-14,'Drizzle coefficient 14 in Y'],
        ['DRZ2Y15', -1.4559746e-13,'Drizzle coefficient 15 in Y'],
        ]

    wcs_keys['direct'] = None
    wcs_keys['dimension'] = [4096, 2048]

    return wcs_keys

def get_WFC3_IR_G102_WCS():
    """
    Defines parameters for the WFC3/IR/G102 slitless mode

    @return: slitless mode parameters
    @rtype: dictionary
    """
    wcs_keys = {}
    #WCS from WCS/WFC3 grism image "ibbu01a3q":
    #/ World Coordinate System and Related Parameters
    wcs_keys['grism'] = [
        ['WCSAXES',                    2, 'number of World Coordinate System axes'],
        ['CRPIX1',   4.920000000000E+02, 'x-coordinate of reference pixel'],
        ['CRPIX2',   5.570000000000E+02, 'y-coordinate of reference pixel'],
        ['CRVAL1',   2.910870382450E+02, 'first axis value at reference pixel'],
        ['CRVAL2',   9.898314109467E+00, 'second axis value at reference pixel'],
        ['CTYPE1', 'RA---TAN'          , 'the coordinate type for the first axis'],
        ['CTYPE2', 'DEC--TAN'          , 'the coordinate type for the second axis'],
        ['CD1_1',         -3.73264E-05, 'partial of first axis coordinate w.r.t. x'],
        ['CD1_2',         -4.60352E-06, 'partial of first axis coordinate w.r.t. y'],
        ['CD2_1',         -5.16424E-06, 'partial of second axis coordinate w.r.t. x'],
        ['CD2_2',          3.34175E-05, 'partial of second axis coordinate w.r.t. y'],
        ['LTV1',        0.0000000E+00, 'offset in X to subsection start'],
        ['LTV2',        0.0000000E+00, 'offset in Y to subsection start'],
        ['LTM1_1',                  1.0, 'reciprocal of sampling rate in X'],
        ['LTM2_2',                  1.0, 'reciprocal of sampling rate in Y'],
        ['PA_APER',             -7.84357, 'Position Angle of reference aperture center (de'],
        ['VAFACTOR',                  1.0, 'velocity aberration plate scale factor'],
        ['ORIENTAT',             -7.84357, 'position angle of image y axis (deg. e of n)'],
        ['RA_APER',   2.910870382450E+02, 'RA of aperture reference position'],
        ['DEC_APER',   9.898314109467E+00, 'Declination of aperture reference position'],
        ['EXPNAME', 'ibbu01a3q' ,'exposure identifier'],
        ]

    wcs_keys['drizzle'] = [
        ['DRZCNUM',15,'Number of coefficients per coordinate'],
        ['DRZSCALE',0.128254,'Scale for drizzling'],
        ['DRZ2X01',           0.14763942, 'Drizzle coefficient 01 in X'],
        ['DRZ2X02',            1.0560001, 'Drizzle coefficient 02 in X'],
        ['DRZ2X03',       -0.00016301818, 'Drizzle coefficient 03 in X'],
        ['DRZ2X04',       -1.9466206E-07, 'Drizzle coefficient 04 in X'],
        ['DRZ2X05',        2.5865694E-05, 'Drizzle coefficient 05 in X'],
        ['DRZ2X06',        7.3464039E-08, 'Drizzle coefficient 06 in X'],
        ['DRZ2X07',       -1.7498839E-10, 'Drizzle coefficient 07 in X'],
        ['DRZ2X08',        1.1485823E-10, 'Drizzle coefficient 08 in X'],
        ['DRZ2X09',        8.1925668E-11, 'Drizzle coefficient 09 in X'],
        ['DRZ2X10',        3.7606885E-11, 'Drizzle coefficient 10 in X'],
        ['DRZ2X11',       -5.1018929E-13, 'Drizzle coefficient 11 in X'],
        ['DRZ2X12',        7.4936782E-13, 'Drizzle coefficient 12 in X'],
        ['DRZ2X13',        3.9890904E-14, 'Drizzle coefficient 13 in X'],
        ['DRZ2X14',        6.9875665E-13, 'Drizzle coefficient 14 in X'],
        ['DRZ2X15',       -3.4250334E-13, 'Drizzle coefficient 15 in X'],
        ['DRZ2Y01',   -9.001424310000001, 'Drizzle coefficient 01 in Y'],
        ['DRZ2Y02',         0.0033880965, 'Drizzle coefficient 02 in Y'],
        ['DRZ2Y03',           0.94300664, 'Drizzle coefficient 03 in Y'],
        ['DRZ2Y04',        6.6289689E-06, 'Drizzle coefficient 04 in Y'],
        ['DRZ2Y05',       -1.2124036E-07, 'Drizzle coefficient 05 in Y'],
        ['DRZ2Y06',        2.8370317E-05, 'Drizzle coefficient 06 in Y'],
        ['DRZ2Y07',        1.7339207E-11, 'Drizzle coefficient 07 in Y'],
        ['DRZ2Y08',       -2.2162965E-10, 'Drizzle coefficient 08 in Y'],
        ['DRZ2Y09',        9.7184953E-12, 'Drizzle coefficient 09 in Y'],
        ['DRZ2Y10',       -1.3682258E-10, 'Drizzle coefficient 10 in Y'],
        ['DRZ2Y11',       -5.4319899E-13, 'Drizzle coefficient 11 in Y'],
        ['DRZ2Y12',        2.3592402E-13, 'Drizzle coefficient 12 in Y'],
        ['DRZ2Y13',        1.4153352E-13, 'Drizzle coefficient 13 in Y'],
        ['DRZ2Y14',       -2.4528943E-14, 'Drizzle coefficient 14 in Y'],
        ['DRZ2Y15',        7.3598547E-13, 'Drizzle coefficient 15 in Y'],
        ]

    wcs_keys['direct'] = None
    wcs_keys['dimension'] = [1014, 1014]

    return wcs_keys

def get_WFC3_UV_G280_WCS():
    """
    Defines parameters for the WFC3/UVI/G280 slitless mode

    @return: slitless mode parameters
    @rtype: dictionary
    """
    wcs_keys = {}
    #WCS from WCS/WFC3 grism image "j8qq12kgq":
    #/ World Coordinate System and Related Parameters
    wcs_keys['grism'] = [
        ['WCSAXES', 2, 'number of World Coordinate System axes'],
        ['CRPIX1', 2.048000000000E+03, 'x-coordinate of reference pixel'],
        ['CRPIX2', -3.920000000000E+02, 'y-coordinate of reference pixel'],
        ['CRVAL1', 1.337668583459E+02, 'first axis value at reference pixel'],
        ['CRVAL2', -4.759016919106E+01, 'second axis value at reference pixel'],
        ['CTYPE1', 'RA---TAN', 'the coordinate type for the first axis'],
        ['CTYPE2', 'DEC--TAN', 'the coordinate type for the second axis'],
        ['CD1_1', -1.02089E-06, 'partial of first axis coordinate w.r.t. x'],
        ['CD1_2', 1.08208E-05, 'partial of first axis coordinate w.r.t. y'],
        ['CD2_1', 1.09748E-05, 'partial of second axis coordinate w.r.t. x'],
        ['CD2_2', 1.76909E-06, 'partial of second axis coordinate w.r.t. y'],
        ['LTV1', 0.0000000E+00, 'offset in X to subsection start'],
        ['LTV2', -6.4200000E+02, 'offset in Y to subsection start'],
        ['LTM1_1', 1.0, 'reciprocal of sampling rate in X'],
        ['LTM2_2', 1.0, 'reciprocal of sampling rate in Y'],
        ['PA_APER', 80.7149, 'Position Angle of reference aperture center (de'],
        ['VAFACTOR', 1.0, 'velocity aberration plate scale factor'],
        ['ORIENTAT', 80.7149, 'position angle of image y axis (deg. e of n)'],
        ['RA_APER', 1.337668583459E+02, 'RA of aperture reference position'],
        ['DEC_APER',  -4.759016919106E+01, 'Declination of aperture reference position'],
        ['EXPNAME','ibbr02bzq' ,'exposure identifier'],
        ]

    wcs_keys['drizzle'] = [
        ['DRZCNUM',15,'Number of coefficients per coordinate'],
        ['DRZSCALE',0.039622,'Scale for drizzling'],
        ['DRZ2X01', -11.42143024, 'Drizzle coefficient 01 in X'],
        ['DRZ2X02', 1.0037689, 'Drizzle coefficient 02 in X'],
        ['DRZ2X03', 0.0018717218, 'Drizzle coefficient 03 in X'],
        ['DRZ2X04', 2.8830796E-06, 'Drizzle coefficient 04 in X'],
        ['DRZ2X05', -2.9770541E-06, 'Drizzle coefficient 05 in X'],
        ['DRZ2X06', 8.511265800000001E-08, 'Drizzle coefficient 06 in X'],
        ['DRZ2X07', 2.0556244E-11, 'Drizzle coefficient 07 in X'],
        ['DRZ2X08', -1.088005E-11, 'Drizzle coefficient 08 in X'],
        ['DRZ2X09', 1.4832995E-11, 'Drizzle coefficient 09 in X'],
        ['DRZ2X10', 2.2900636E-11, 'Drizzle coefficient 10 in X'],
        ['DRZ2X11', 1.6923456E-15, 'Drizzle coefficient 11 in X'],
        ['DRZ2X12', 7.120595E-16, 'Drizzle coefficient 12 in X'],
        ['DRZ2X13', -1.7063885E-14, 'Drizzle coefficient 13 in X'],
        ['DRZ2X14', -4.3558507E-15, 'Drizzle coefficient 14 in X'],
        ['DRZ2X15', -1.5479474E-14, 'Drizzle coefficient 15 in X'],
        ['DRZ2Y01', -2.23386146, 'Drizzle coefficient 01 in Y'],
        ['DRZ2Y02', 0.061529103, 'Drizzle coefficient 02 in Y'],
        ['DRZ2Y03', 1.0054915, 'Drizzle coefficient 03 in Y'],
        ['DRZ2Y04', 1.3810195E-07, 'Drizzle coefficient 04 in Y'],
        ['DRZ2Y05', 2.6551435E-06, 'Drizzle coefficient 05 in Y'],
        ['DRZ2Y06', -3.0875627E-06, 'Drizzle coefficient 06 in Y'],
        ['DRZ2Y07', 3.6994702E-12, 'Drizzle coefficient 07 in Y'],
        ['DRZ2Y08', 1.62375E-11, 'Drizzle coefficient 08 in Y'],
        ['DRZ2Y09', -1.0216096E-11, 'Drizzle coefficient 09 in Y'],
        ['DRZ2Y10', 1.0595473E-11, 'Drizzle coefficient 10 in Y'],
        ['DRZ2Y11', 6.5447373E-16, 'Drizzle coefficient 11 in Y'],
        ['DRZ2Y12', 1.2674826E-15, 'Drizzle coefficient 12 in Y'],
        ['DRZ2Y13', 1.1562322E-14, 'Drizzle coefficient 13 in Y'],
        ['DRZ2Y14', -8.6793139E-15, 'Drizzle coefficient 14 in Y'],
        ['DRZ2Y15', -1.3467473E-15, 'Drizzle coefficient 15 in Y'],
        ]

    wcs_keys['direct'] = None
    wcs_keys['dimension'] = [4096, 2048]

    return wcs_keys

def get_NICMOS3_G141_WCS():
    """
    Defines parameters for the NICMOS/G141 slitless mode

    @return: slitless mode parameters
    @rtype: dictionary
    """
    wcs_keys = {}

    # WCS from NICMOS G141 image n6le01upq
    #/ World Coordinate System and Related Parameters
    wcs_keys['grism'] = [
        ['WCSAXES',2,'number of World Coordinate System axes'],
        ['CRPIX1',128.0,'x-coordinate of reference pixel'],
        ['CRPIX2',128.0,'y-coordinate of reference pixel'],
        ['CRVAL1',1.291877276104E+02,'first axis value at reference pixel'],
        ['CRVAL2',9.136572816451E-01,'second axis value at reference pixel'],
        ['CTYPE1','RA---TAN','the coordinate type for the first axis'],
        ['CTYPE2','DEC--TAN','the coordinate type for the second axis'],
        ['CD1_1',-3.01591E-05,'partial of first axis coordinate w.r.t. x'],
        ['CD1_2',4.7564E-05,'partial of first axis coordinate w.r.t. y'],
        ['CD2_1',4.76877E-05,'partial of second axis coordinate w.r.t. x'],
        ['CD2_2',3.00809E-05,'partial of second axis coordinate w.r.t. y'],
        ]

    wcs_keys['drizzle'] = [
        ['DRZCNUM',                  10,'Number of coefficients per coordinate'],
        ['DRZSCALE', 0.2,'Scale for drizzling'],
        ['DRZ2X01',  0.0,'Drizzle coefficient 01 in X'],
        ['DRZ2X02',  1.0018288,'Drizzle coefficient 02 in X'],
        ['DRZ2X03',  0.0,'Drizzle coefficient 03 in X'],
        ['DRZ2X04',  8.034670000000001E-06,'Drizzle coefficient 04 in X'],
        ['DRZ2X05',  1.32241E-05,'Drizzle coefficient 05 in X'],
        ['DRZ2X06',  5.83064E-06,'Drizzle coefficient 06 in X'],
        ['DRZ2X07',  0.0,'Drizzle coefficient 07 in X'],
        ['DRZ2X08',  0.0,'Drizzle coefficient 08 in X'],
        ['DRZ2X09',  0.0,'Drizzle coefficient 09 in X'],
        ['DRZ2X10',  0.0,'Drizzle coefficient 10 in X'],
        ['DRZ2Y01',  0.0,'Drizzle coefficient 01 in Y'],
        ['DRZ2Y02', -0.000893359,'Drizzle coefficient 02 in Y'],
        ['DRZ2Y03',  0.99816635,'Drizzle coefficient 03 in Y'],
        ['DRZ2Y04', -1.80668E-05,'Drizzle coefficient 04 in Y'],
        ['DRZ2Y05',  5.989E-07,'Drizzle coefficient 05 in Y'],
        ['DRZ2Y06', -1.15787E-05,'Drizzle coefficient 06 in Y'],
        ['DRZ2Y07', 0.0,'Drizzle coefficient 07 in Y'],
        ['DRZ2Y08', 0.0,'Drizzle coefficient 08 in Y'],
        ['DRZ2Y09', 0.0,'Drizzle coefficient 09 in Y'],
        ['DRZ2Y10', 0.0,'Drizzle coefficient 10 in Y']
        ]

    wcs_keys['direct'] = None
    wcs_keys['dimension'] = [256, 256]

    return wcs_keys

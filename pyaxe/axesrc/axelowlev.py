from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import os
import subprocess

from .axeerror import aXeError
from . import axeutils

# define the good return
# value for the binaries
GOOD_RETURN_VALUE = 0


class TaskWrapper(object):
    """General class to execute C-tasks"""
    def __init__(self, taskname="", tshort=""):
        """

        Parameters
        ----------
        taskname: str
            name of the C-executable
        tshort: str
            shor name for task
        """
        self.taskname = taskname
        self.tshort = tshort

        # initialize the command list
        self.command_list = []

        # save a name for stdout
        self.stdout = axeutils.getOUTPUT(tshort+'.stdout')

        # save a name for stderr
        self.stderr = axeutils.getOUTPUT(tshort+'.stderr')

        # put the command into the list
        self.command_list.append(axeutils.getBINDIR(taskname))

    def _cleanup(self):
        """The method deletes the files created for stdout and stderr.

        This is a usual cleaning procedure in case nothing bad happened.
        """
        # delete stdout/stderr
        if os.path.isfile(self.stdout):
            os.unlink(self.stdout)
        if os.path.isfile(self.stderr):
            os.unlink(self.stderr)

    def _report_all(self, silent=True):
        """Print stdout and stderr on the screen

        The method gives a feedback in case of problems. stdout and stderr
        are both listed onto the screen for a further interactive analysis.

        Parameters
        ----------
        silent: bool
            indicates silent/noisy runs
        """
        # check whether the command
        # was run silent
        if silent:
            # dump the files with
            # stdout and stderr onto the screen
            print("\nThere was a problem in the task: {0:s}"
                  .format(self.taskname))
            print("The output of the task ({0:s}) is:\n".format(self.stdout))
            print("---------------------------------------------------------"
                  "-----------------------")
            for line in open(self.stdout):
                print(line.strip())
            print("\n\nThe error report is of the task ({0:s}) is\n:"
                  .format(self.stdout))
            print("----------------------------------------------------------"
                  "----------------------")
            for line in open(self.stderr):
                print(line.strip())

        # report an error
        raise aXeError("An error occurred in the aXe task: {0:s}"
                       .format(self.taskname))

    def run(self, silent=False):
        """Run the wrapped task

        The method executes the associated C-executable. The return code given
        by the C-executable is returned. In silent mode stdout and stderr
        are writtren to a file, in non-silent mode to the screen.

        Parameters
        ----------
        silent: bool
            for silent mode

        Returns
        -------
        retcode: int
            the return code of the C-executable
        """
        # is output desired
        if silent:
            # open stdout/stderr
            sout = open(self.stdout, 'w+')
            serr = open(self.stderr, 'w+')

            # execute the task
            retcode = subprocess.call(self.command_list,
                                      stdout=sout,
                                      stderr=serr)

            # close stdout/stderr
            sout.close()
            serr.close()

        else:

            # execute the task with the default stdout and
            # stderr, which is the system one
            print(self.command_list)
            retcode = subprocess.call(self.command_list)

        # return the result
        return retcode

    def runall(self, silent=False):
        """Run the wrapped task

        The method executes the associated C-executable. The return code given
        by the C-executable is returned. In silent mode stdout and stderr
        are writtren to a file, in non-silent mode to the screen.

        Parameters
        ----------
        silent: bool
            for silent mode

        Returns
        -------
        recode: int
            he return code of the C-executable
        """
        # run the executable
        retcode = self.run(silent=silent)

        # check whether the run was good
        if retcode == GOOD_RETURN_VALUE:
            # do the cleaning
            self._cleanup()
        else:
            self._report_all(silent)

        # return the result
        return retcode


class aXe_AF2PET(TaskWrapper):
    """Wrapper around the aXe_AF2PET task"""
    def __init__(self, grism, config, **params):
        """This method is a simple initializer for the class.

        All variables are transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_AF2PET, self).__init__('aXe_AF2PET', 'af2pet')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for background flag
        if (('back' in params) and (params['back'])):
            # put the bck-flag to the list
            self.command_list.append('-bck')

        # check for the in-AF name
        if (('in_af' in params) and (params['in_af'] is not None)):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # check for the out-PET name
        if (('out_pet' in params) and (params['out_pet'] is not None)):
            # put the name to the list
            self.command_list.append("-out_PET={0:s}"
                                     .format(params['out_pet']))


class aXe_GPS(TaskWrapper):
    """Wrapper around the aXe_GPS task"""
    def __init__(self, grism, config, beam_ref, xval, yval):
        """This method is a simple initializer for the class.

        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        beam_ref: str
            beam to use
        xval: float
            x-position to check
        yval: float
            y-position to check
        """
        # initialize via superclass
        super(aXe_GPS, self).__init__('aXe_GPS', 'axegps')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # put the beam name to the list
        self.command_list.append(beam_ref)

        # put the x-value to the list
        self.command_list.append(str(xval))

        # put the y-value to the list
        self.command_list.append(str(yval))


class aXe_BE(TaskWrapper):
    """Wrapper around the aXe_BE task"""
    def __init__(self, grism, config, **params):
        """This method is a simple initializer for the class.

        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_BE, self).__init__('aXe_BE', 'backest')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'np'
        if (('np' in params) and (params['np'] is not None)):
            self.command_list.append('-np={0:s}'.format(str(params['np'])))

        # append the parameter 'interp'
        if (('interp' in params) and (params['interp'] is not None)):
            self.command_list.append("-interp={0:s}"
                                     .format(str(params['interp'])))

        # append the parameter 'niter_med'
        if (('niter_med' in params) and (params['niter_med'] is not None)):
            self.command_list.append('-niter_med={0:s}'
                                     .format(str(params['niter_med'])))

        # append the parameter 'niter_fit'
        if (('niter_fit' in params) and (params['niter_fit'] is not None)):
            self.command_list.append('-niter_fit={0:s}'
                                     .format(str(params['niter_fit'])))

        # append the parameter 'kappa'
        if (('kappa' in params) and (params['kappa'] is not None)):
            self.command_list.append('-kappa={0:s}'
                                     .format(str(params['kappa'])))

        # append the parameter 'smooth_length'
        if (params['smooth_length'] is not None):
            self.command_list.append('-smooth_length={0:s}'
                                     .format(str(params['smooth_length'])))

        # append the parameter 'smooth_length'
        if (params['smooth_fwhm'] is not None):
            self.command_list.append('-fwhm={0:s}'
                                     .format(str(params['smooth_fwhm'])))

        # append the parameter 'af_file'
        if (params['in_af'] is not None):
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # append the parameter 'out_bck'
        if (params['out_bck'] is not None):
            self.command_list.append('-out_BCK={0:s}'
                                     .format(params['out_bck']))

        # append the flag 'old_bck'
        if (('old_bck' in params) and (params['old_bck'])):
            self.command_list.append('-nor_flag')

        # append the flag 'mask'
        if (('mask' in params) and (params['mask'])):
            self.command_list.append('-msk')


class aXe_DRZ2PET(TaskWrapper):
    """Wrapper around the aXe_DRZ2PET task"""
    def __init__(self, inlist, config, **params):
        """This method is a simple initializer for the class.

        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_DRZ2PET, self).__init__('aXe_DRZ2PET', 'drz2pet')

        # put the grism name to the list
        self.command_list.append(inlist)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'out_pet'
        if (('out_pet' in params) and (params['out_pet'] is not None)):
            # put the name to the list
            self.command_list.append('-out_PET={0:s}'
                                     .format(params['out_pet']))

        # append the parameter 'in_af'
        if (('in_af' in params) and (params['in_af'] is not None)):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # append the flag 'bck'
        if (('back' in params) and (params['back'])):
            # put the bck-flag to the list
            self.command_list.append('-bck')

        # append the flag 'bck'
        if (('opt_extr' in params) and (params['opt_extr'])):
            # put the bck-flag to the list
            self.command_list.append('-opt_extr')


class aXe_DRZPREP(TaskWrapper):
    """Wrapper around the aXe_DRZPREP task"""
    def __init__(self, inlist, configs, **params):
        """This method is a simple initializer for the class.
        All variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        inlist: str
            name of the input image list
        configs: str
            aXe configuration file term
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_DRZPREP, self).__init__('aXe_DRZPREP', 'drzprep')

        # store the 'back'-flag
        if (('back' in params) and (params['back'])):
            # put the opt_extr-flag to the list
            self.bck = True
        else:
            self.bck = False

        # put the grism name to the list
        self.command_list.append(inlist)

        # put the config file name to the list
        self.command_list.append(configs)

        # append the flag 'opt_extr'
        if (('opt_extr' in params) and (params['opt_extr'])):
            # put the opt_extr-flag to the list
            self.command_list.append('-opt_extr')

    def runall(self, silent=False):
        """Run the wrapped task

        The method executes the associated C-executable. The return code given
        by the C-executable is returned. In silent mode stdout and stderr
        are writtren to a file, in non-silent mode to the screen.

        Parameters
        ----------
        silent: bool
            or silent mode

        Returns
        -------
        recode: int
            the return code of the C-executable
        """
        # run the method of the super class
        super(aXe_DRZPREP, self).runall(silent=silent)

        # check for the background flag
        if self.bck:
            # put the bck-flag to the list
            self.command_list.append('-bck')

            # run the method of the super class
            super(aXe_DRZPREP, self).runall(silent=silent)


class aXe_GOL2AF(TaskWrapper):
    """Wrapper around the aXe_GOL2AF task"""
    def __init__(self, grism, config, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        configs: str
            aXe configuration file term
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_GOL2AF, self).__init__('aXe_GOL2AF', 'gol2af')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'mfwhm'
        if (params['mfwhm'] is not None):
            # put the name to the list
            self.command_list.append('-mfwhm={0:s}'
                                     .format(str(params['mfwhm'])))

        # append the parameter 'dmag'
        if (params['dmag'] is not None):
            # put the name to the list
            self.command_list.append('-dmag={0:s}'.format(str(params['dmag'])))

        # append the parameter 'lambda_mark'
        if (params['lambda_mark'] is not None):
            # put the name to the list
            self.command_list.append('-lambda_mark={0:s}'
                                     .format(str(params['lambda_mark'])))
        else:
            self.command_list.append('-lambda_mark=800.0')

        # append the parameter 'out_pet'
        if (params['out_pet'] is not None):
            # put the name to the list
            self.command_list.append('-out_PET={0:s}'
                                     .format(params['out_pet']))

        # append the parameter 'in_af'
        if (params['out_af'] is not None):
            # put the name to the list
            self.command_list.append('-out_AF={0:s}'.format(params['out_af']))

        # append the parameter 'in_gol'
        if (params['in_gol'] is not None):
            # put the name to the list
            self.command_list.append('-in_GOL={0:s}'.format(params['in_gol']))

        # append the flag 'slitless_geom'
        if (params['slitless_geom']):
            # put the exclude_faint-flag to the list
            self.command_list.append('-slitless_geom=1')
        else:
            self.command_list.append('-slitless_geom=0')

        # append the flag 'orient'
        if (params['orient']):
            # put the exclude_faint-flag to the list
            self.command_list.append('-orient=1')
        else:
            self.command_list.append('-orient=0')

        # append the flag 'exclude'
        if (params['exclude']):
            # put the exclude_faint-flag to the list
            self.command_list.append('-exclude_faint')

        #  append the flag 'bck'
        if (params['back']):
            # put the bck-flag to the list
            self.command_list.append('-bck')


class aXe_INTPIXCORR(TaskWrapper):
    """Wrapper around the aXe_INTPIXCORR task"""
    def __init__(self, grism, config, **params):
        """
        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_INTPIXCORR, self).__init__('aXe_INTPIXCORR', 'ipixcorr')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_OAF' name
        if (params['in_OAF'] is not None):
            # put the name to the list
            self.command_list.append('-in_OAF={0:s}'
                                     .format(str(params['in_oaf'])))

        # check for the 'in_SPC' name
        if (params['in_SPC'] is not None):
            # put the name to the list
            self.command_list.append('-in_SPC={0:s}'
                                     .format(str(params['in_spc'])))

        # check for the 'out_SPC' name
        if (params['out_SPC'] is not None):
            # put the name to the list
            self.command_list.append('-out_SPC={0:s}'
                                     .format(params['out_spc']))

        # check for the 'max_ext' value
        if (params['max_ext'] is not None):
            # put the value to the list
            self.command_list.append('-max_ext={0:s}'
                                     .format(str(params['max_ext'])))
        else:
            # put the default to the list
            self.command_list.append('-max_ext=0.0')

        # append the flag 'ip_corr'
        if (params['ip_corr']):
            # put the ip_corr-flag to the list
            self.command_list.append('-ipixcorr')

        # append the flag 'nl_corr'
        if (params['nl_corr']):
            # put the nl_corr-flag to the list
            self.command_list.append('-nlincorr')


class aXe_NICBACK(TaskWrapper):
    """Wrapper"""
    def __init__(self, grism, config, master_bck, back_ped=None):
        """
        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_NICBACK, self).__init__('aXe_NICBACK', 'nicback')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # put the master background
        self.command_list.append(master_bck)

        if ((back_ped is not None) and (len(back_ped) > 1)):
            # put the master background
            self.command_list.append(back_ped)


class aXe_SCALEBCK(TaskWrapper):
    """
    Wrapper around the aXe_PET2SPC task
    """
    def __init__(self, grism, mask, config, master_sky, to_master=False,
                 make_plis=False):
        """
        Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        @param grism: name of the dispersed image
        @type grism: string
        @param mask: name of the mask image
        @type mask: string
        @param config: name of the aXe configuration file
        @type config: string
        @param master_sky: name of the master sky image
        @type master_sky: string
        @param to_master: scale grism to master
        @type to_master: boolean
        @param make_plis: generate the pixel list
        @type make_plis: boolean
        """
        # initialize via superclass
        super(aXe_SCALEBCK, self).__init__('aXe_SCALEBCK', 'scalebck')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the mask name to the list
        self.command_list.append(mask)

        # put the config file name to the list
        self.command_list.append(config)

        # put the master background
        self.command_list.append(master_sky)

        if to_master:
            # add the flag for scaling to the master image
            self.command_list.append('-toMaster')

        if make_plis:
            # add flag to generate pixel list
            self.command_list.append('-make_plis')


class aXe_PET2SPC(TaskWrapper):
    """Wrapper around the aXe_PET2SPC task"""
    def __init__(self, grism, config, **params):
        """
        Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_PET2SPC, self).__init__('aXe_PET2SPC', 'pet2spc')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_AF' name
        if (('in_af' in params) and (params['in_af'] is not None)):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # check for the 'OPET' name
        if (params['opet'] is not None):
            # put the name to the list
            self.command_list.append('-OPET={0:s}'.format(params['opet']))

        # check for the 'BPET' name
        if (params['bpet'] is not None):
            # put the name to the list
            self.command_list.append('-BPET={0:s}'.format(params['bpet']))

        # check for the 'out_SPC' name
        if (params['out_spc'] is not None):
            # put the name to the list
            self.command_list.append('-out_SPC={0:s}'
                                     .format(params['out_spc']))

        # append the flag 'drzpath'
        if (params['drzpath']):
            # put the ip_corr-flag to the list
            self.command_list.append('-drz')

        # append the flag 'noBPET'
        if (not params['use_bpet']):
            # put the ip_corr-flag to the list
            self.command_list.append('-noBPET')

        # append the flag 'opt_weights'
        if (params['weights']):
            # put the ip_corr-flag to the list
            self.command_list.append('-opt_weights')

        # append the flag 'noflux'
        if (not params['do_flux']):
            # put the ip_corr-flag to the list
            self.command_list.append('-noflux')

        # append the flag 'smooth_conv'
        if (('adj_sens' in params) and (params['adj_sens'])):
            # put the ip_corr-flag to the list
            self.command_list.append('-smooth_conv')


class aXe_PETCONT(TaskWrapper):
    """Wrapper around the aXe_PETCONT task"""
    def __init__(self, grism, config, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters        """
        # initialize via superclass
        super(aXe_PETCONT, self).__init__('aXe_PETCONT', 'petcont')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_AF' name
        if (('in_af' in params) and (params['in_af'] is not None)):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # check for the 'specmodels' name
        if (params['spec_models'] is not None):
            # put the name to the list
            self.command_list.append('-model_spectra={0:s}'
                                     .format(params['spec_models']))

        # check for the 'objectmodels' name
        if (params['object_models'] is not None):
            # put the name to the list
            self.command_list.append('-model_images={0:s}'
                                     .format(params['object_models']))

        # append the flag 'cont_model'
        if (params['cont_model'] is 'gauss'):
            # put the gauss-flag to the list
            self.command_list.append('-cont_model=1')
        elif (params['cont_model'] is 'direct'):
            # put the according number to the list
            self.command_list.append('-cont_model=2')
        elif (params['cont_model'] is 'fluxcube'):
            # put the according number to the list
            self.command_list.append('-cont_model=3')
        elif (params['cont_model'] is 'geometric'):
            # put the according number to the list
            self.command_list.append('-cont_model=4')

        # append the flag 'inter_type'
        if (params['inter_type'] is 'linear'):
            # put the according number to the list
            self.command_list.append('-inter_type=1')
        elif (params['inter_type'] is 'polynomial'):
            # put the according number to the list
            self.command_list.append('-inter_type=2')
        elif (params['inter_type'] is 'spline'):
            # put the according number to the list
            self.command_list.append('-inter_type=3')

        # append the parameter 'model_scale'
        if (params['model_scale'] is not None):
            # put the flag to the list
            self.command_list.append('-model_scale={0:s}'
                                     .format(str(params['model_scale'])))

        # append the parameter 'lambda_psf'
        if (params['lambda_psf'] is not None):
            # put the number to the list
            self.command_list.append('-lambda_psf={0:s}'
                                     .format(str(params['lambda_psf'])))
        else:
            self.command_list.append('-lambda_psf=800.0')

        # append the flag 'cont_map'
        if (params['cont_map']):
            # put the cont_map-flag to the list
            self.command_list.append('-cont_map')

        # append the flag 'no_pet' to indicate
        # that no PET exists
        if (params['no_pet']):
            # append the no-PET flagg
            self.command_list.append('-noPET')


class aXe_PETFF(TaskWrapper):
    """Wrapper around the aXe_PETFF task"""
    def __init__(self, grism, config, **params):
        """
        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_PETFF, self).__init__('aXe_PETFF', 'petff')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'FFNAME' name
        if (params['ffname'] is not None):
            # put the name to the list
            self.command_list.append('-FFNAME={0:s}'.format(params['ffname']))

        # append the flag 'bck'
        if (params['back']):
            # put the according flag to the list
            self.command_list.append('-bck')


class aXe_PETIPC(TaskWrapper):
    """Wrapper around the aXe_PETIPC task"""
    def __init__(self, grism, config, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_PETIPC, self).__init__('aXe_PETIPC', 'petipc')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # check for the 'in_OAF' name
        if (('in_OAF' in params) and (params['in_OAF'] is not None)):
            # put the name to the list
            self.command_list.append('-in_OAF={0:s}'.format(params['in_oaf']))

        # check for the 'in_PET' name
        if (('in_PET' in params) and (params['in_PET'] is not None)):
            # put the name to the list
            self.command_list.append('-in_PET={0:s}'.format(params['in_pet']))

        # check for the 'out_PET' name
        if (('out_PET' in params) and (params['out_PET'] is not None)):
            # put the name to the list
            self.command_list.append('-out_PET={0:s}'
                                     .format(params['out_pet']))

        # append the parameter 'max_ext'
        if (('max_ext' in params) and (params['max_ext'] is not None)):
            # put the number to the list
            self.command_list.append('-max_ext={0:s}'
                                     .format(str(params['max_ext'])))

        # append the flag 'bck'
        if (params['back']):
            # put the according flag to the list
            self.command_list.append('-bck')

        # append the flag 'origname'
        if (('orig_name' in params) and (params['orig_name'])):
            # put the according flag to the list
            self.command_list.append('-origname')


class aXe_SEX2GOL(TaskWrapper):
    """Wrapper around the aXe_SEX2GOL task"""
    def __init__(self, grism, config, in_sex, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        iolname: str
            name of the input object list
        dirname: str
            name of the direct image

        """
        # initialize via superclass
        super(aXe_SEX2GOL, self).__init__('aXe_SEX2GOL', 'sex2gol')

        # check whether a direct image exists
        if (('use_direct' in params) and (params['use_direct'])):
            # put the direct image name to the list
            self.command_list.append(params['dirname'])

            # put the grism name to the list
            self.command_list.append(grism)

            # put the grism name to the list
            self.command_list.append(config)
        else:
            # put the grism name to the list
            self.command_list.append(grism)

            # put the grism name to the list
            self.command_list.append(config)

            # mark that there is no direct image
            self.command_list.append('-no_direct_image')

        # put the SExtractor cat to the list
        self.command_list.append('-in_SEX={0:s}' % in_sex)

        # check for the 'out_SEX' name
        if (params['out_sex'] is not None):
            # put the name to the list
            self.command_list.append('-out_SEX={0:s}'
                                     .format(params['out_sex']))

        # append the parameter 'dir_hdu'
        if (params['dir_hdu'] is not None):
            # put the number to the list
            self.command_list.append('-dir_hdu={0:s}'
                                     .format(str(params['dir_hdu'])))

        # append the parameter 'spec_hdu'
        if (params['spec_hdu'] is not None):
            # put the number to the list
            self.command_list.append('-spec_hdu={0:s}'
                                     .format(str(params['spec_hdu'])))


class aXe_STAMPS(TaskWrapper):
    """Wrapper around the aXe_STAMPS task"""
    def __init__(self, grism, config, **params):
        """
        Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        @param grism: name of the dispersed image
        @type grism: string
        @param config: name of the aXe configuration file
        @type config: string
        @param **params: all other parameters
        @type **params: dictionary
        """
        # initialize via superclass
        super(aXe_STAMPS, self).__init__('aXe_STAMPS', 'stamps')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'in_af'
        if (params['in_af'] is not None):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))

        # check for the 'in_PET' name
        if (params['in_pet'] is not None):
            # put the name to the list
            self.command_list.append('-in_PET={0:s}'.format(params['in_pet']))

        # check for the 'out_STP' name
        if (params['out_stp'] is not None):
            # put the name to the list
            self.command_list.append('-out_STP={0:s}'
                                     .format(params['out_stp']))

        # append the flag 'drz'
        if (params['sampling'] is 'drizzle'):
            # put the flagg to the list
            self.command_list.append('-drzstamp')
        elif (params['sampling'] is 'rectified'):
            # put the flagg to the list
            self.command_list.append('-rectified')

        # append the flag 'drz'
        if (params['drzpath']):
            # put the according flag to the list
            self.command_list.append('-drz')


class aXe_TRACEFIT(TaskWrapper):
    """Wrapper around the aXe_TRACEFIT task"""
    def __init__(self, grism, config, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        grism: str
            name of the dispersed image
        config: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_TRACEFIT, self).__init__('aXe_TRACEFIT', 'tracefit')

        # put the grism name to the list
        self.command_list.append(grism)

        # put the config file name to the list
        self.command_list.append(config)

        # append the parameter 'in_af'
        if (('in_af' in params) and (params['in_af'] is not None)):
            # put the name to the list
            self.command_list.append('-in_AF={0:s}'.format(params['in_af']))


class aXe_FILET(TaskWrapper):
    """Wrapper around the aXe_STAMPS task"""
    def __init__(self, dppfile, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        @param dppfile: name of the dpp-file
        @type dppfile: string
        @param **params: all other parameters
        @type **params: dictionary
        """
        # initialize via superclass
        super(aXe_FILET, self).__init__('aXe_FILET', 'filet')

        # put the grism name to the list
        self.command_list.append(dppfile)

        # append the 'opt_extr' flag
        if (('opt_extr' in params) and (params['opt_extr'])):
            self.command_list.append('-opt_extr')

            # append the parameter 'in_af'
        if (('drztmp' in params) and (params['drztmp'] is not None)):
            # put the name to the list
            self.command_list.append('-drztmp={0:s}'.format(params['drztmp']))


class aXe_DIRIMAGE(TaskWrapper):
    """Wrapper around the aXe_DIRIMAGE task"""
    def __init__(self, dirname, config, tpass_direct, **params):
        """Initializer for the class

        This method is a simple initializer for the class. All
        variables a transferred to a list, if necessary with the
        appropriate leading parameter name

        Parameters
        ----------
        dirname: str
            name of the direct image
        configfile: str
            name of the aXe configuration file
        params: dict
            all other parameters
        """
        # initialize via superclass
        super(aXe_DIRIMAGE, self).__init__('aXe_DIRIMAGE', 'dirimage')

        # put the direct image name to the list
        self.command_list.append(dirname)

        # put the aXe configuration file name to the list
        self.command_list.append(config)

        # put the total passband file name to the list
        self.command_list.append(tpass_direct)

        # check whether model spectra are given
        # append the file name to the list
        if (params['model_spectra'] is not None):
            # self.command_list.append('-model_spectra='+str(model_spectra))
            self.command_list.append('-model_spectra={0:s}'
                                     .format(params['model_spectra']))

        # check whether model images are given
        # append the file name to the list
        if (params['model_images'] is not None):
            self.command_list.append('-model_images={0:s}'
                                     .format(params['model_images']))

        # check whether model images are given
        # append the file name to the list
        if (params['tel_area'] is not None):
            self.command_list.append('-tel_area={0:f}'
                                     .format(params['tel_area']))

        # append the model scale, give default if necessary
        if (params['model_scale'] is not None):
            self.command_list.append('-model_scale={0:f}'
                                     .format(params['model_scale']))
        else:
            self.command_list.append('-model_scale=5.0')

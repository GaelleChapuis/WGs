import os
from ibllib.io import params as iopar
from getpass import getpass
from pathlib import Path, PurePath

_PAR_ID_STR = 'one_params'


def default():
    par = {"ALYX_LOGIN": "test_user",
           "ALYX_PWD": None,
           "ALYX_URL": "https://test.alyx.internationalbrainlab.org",
           "CACHE_DIR": str(PurePath(Path.home(), "Downloads", "FlatIron")),
           "FTP_DATA_SERVER": "ftp://ibl.flatironinstitute.org",
           "FTP_DATA_SERVER_LOGIN": "iblftp",
           "FTP_DATA_SERVER_PWD": None,
           "HTTP_DATA_SERVER": "http://ibl.flatironinstitute.org",
           "HTTP_DATA_SERVER_LOGIN": "iblmember",
           "HTTP_DATA_SERVER_PWD": None,
           "GLOBUS_CLIENT_ID": None,
           }
    return iopar.from_dict(par)


# first get current and default parameters
def setup():
    par_current = iopar.read(_PAR_ID_STR)
    par_default = default()
    if par_current is None:
        par_current = par_default

    def _get_current_par(k):
        cpar = getattr(par_current, k, None)
        if cpar is None:
            cpar = getattr(par_default, k, None)
        return cpar

    par = iopar.as_dict(par_current)
    for k in par.keys():
        cpar = _get_current_par(k)
        if "PWD" not in k:
            par[k] = input("Param " + k + ",  current value is [" + str(cpar) + "]:") or cpar

    cpar = _get_current_par("ALYX_PWD")
    prompt = "Enter the Alyx password for " + par["ALYX_LOGIN"] + '(leave empty to keep current):'
    par["ALYX_PWD"] = getpass(prompt) or cpar

    cpar = _get_current_par("HTTP_DATA_SERVER_PWD")
    prompt = "Enter the FlatIron HTTP password for " + par["HTTP_DATA_SERVER_LOGIN"] +\
             '(leave empty to keep current): '
    par["HTTP_DATA_SERVER_PWD"] = getpass(prompt) or cpar

    cpar = _get_current_par("FTP_DATA_SERVER_PWD")
    prompt = "Enter the FlatIron FTP password for " + par["FTP_DATA_SERVER_LOGIN"] +\
             '(leave empty to keep current): '
    par["FTP_DATA_SERVER_PWD"] = getpass(prompt) or cpar

    # default to home dir if empty dir somehow made it here
    if len(par['CACHE_DIR']) == 0:
        par['CACHE_DIR'] = str(PurePath(Path.home(), "Downloads", "FlatIron"))

    par = iopar.from_dict(par)

    # create directory if needed
    if par.CACHE_DIR and not os.path.isdir(par.CACHE_DIR):
        os.mkdir(par.CACHE_DIR)
    iopar.write(_PAR_ID_STR, par)
    print('ONE Parameter file location: ' + iopar.getfile(_PAR_ID_STR))


def get():
    par = iopar.read(_PAR_ID_STR)
    if par is None:
        setup()
    return iopar.read(_PAR_ID_STR, default=default())

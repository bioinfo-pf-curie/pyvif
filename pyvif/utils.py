# -*- coding: utf-8 -*-

"""Utils fonctions module."""

import base64
import io
import os
import pkg_resources

from pyvif import logger


def include_file(name, b64=False):
    try:
        if b64:
            with io.open(name, "rb") as f:
                return base64.b64encode(f.read()).decode('utf-8')
        else:
            with io.open(name, "r", encoding='utf-8') as f:
                return f.read()
    except (OSError, IOError) as e:
        logger.error("Could not include file '{}': {}".format(name, e))


def template_data(filename=None):
    """ Return full path of a pyvif template data file.

    :param str filename: a valid filename to be found in `templates`.
    :return: the file path.
    """
    pyvif_path = pkg_resources.get_distribution('pyvif').location
    data_dir = os.sep.join([pyvif_path, 'pyvif', 'templates'])
    filename = os.sep.join([data_dir, filename])
    if os.path.exists(filename):
        return filename
    raise FileNotFoundError("Unknown file {}.".format(filename))


def embed_png(plot_function, output_arg, **kwargs):
    """ Take plot function as input and return image as b64 string.
    You must set output argument of plot function to connect the buffer.

    :param func plot_function: plot function to embed.
    :param str output_arg: output argument of plot_function.
    :param **kwargs: plot functions argument.
    """
    buf = io.BytesIO()
    kwargs = dict({output_arg: buf}, **kwargs)
    plot_function(**kwargs)
    return base64.b64encode(buf.getvalue()).decode('utf-8')

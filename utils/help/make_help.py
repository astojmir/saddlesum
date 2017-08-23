#! /usr/bin/env python
#
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Code author:  Aleksandar Stojmirovic
#

"""
Produces help strings for saddlesum and saddlesum-etd from rst documentation.

SYNOPSIS:

    make_help.py rst_file

"""

import sys
import os
import subprocess
import shutil
from sphinx.application import Sphinx
from qmbpmn.common.utils.filesys import makedirs2

help_files = ['saddlesum-cli', 'saddlesum-show-etd']
macro_names = ['HELP_SADDLESUM', 'HELP_SHOW_ETD']


def generate_docs(srcdir):
    """
    Use Sphinx to generate documentation pages from ReST input.
    """


    scriptdir = os.path.abspath(os.path.split(__file__)[0])
    outdir = os.path.join(scriptdir, '.tmp_help')
    templates_path = scriptdir
    in_paths = [os.path.join(srcdir, '%s.rst' % f) for f in help_files]
    out_paths = [os.path.join(outdir, '%s.html' % f) for f in help_files]

    makedirs2(outdir)
    doctreedir = os.path.join(outdir, 'doctrees')
    app = Sphinx(srcdir, srcdir, outdir, doctreedir, 'html',
                 {'templates_path': [templates_path]},
                 None, None, True)
    app.builder.add_permalinks = False
    app.builder.templates.environment.globals['embedded'] = False
    app.builder.build_specific(in_paths)

    # Pass html through elinks to render text
    cmd = 'elinks -force-html -dump'
    args = cmd.split(' ')

    help_texts = []
    for name, f in zip(macro_names, out_paths):
        _args = args + [f]
        process = subprocess.Popen(_args, shell=False,
                                   stdin=None,
                                   stdout=subprocess.PIPE)
        help_texts.append((name, process.communicate()[0]))

    shutil.rmtree(outdir, True)
    return help_texts


def write_c_file(texts, outfile):

    header = '/* This file was automatically generated from RST docs */ \n\n'
    with open(outfile, 'wb') as fp:
        fp.write(header)
        for name, txt in texts:
            lines = txt.splitlines()[2:]
            fp.write('#define %s "" \\\n' % name)
            for line in lines:
                fp.write('"%s\\n" \\\n' % line)
            fp.write('"\\n" \n\n\n')


if __name__ == '__main__':

    srcdir = os.path.abspath(sys.argv[1])
    outfile = os.path.abspath(sys.argv[2])
    texts = generate_docs(srcdir)
    write_c_file(texts, outfile)

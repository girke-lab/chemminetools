#!/usr/bin/python
# -*- coding: utf-8 -*-

"""iterator over an sdf files"""

from subprocess import Popen, PIPE
import re


class GzipFile(object):

    def __init__(self, path):
        self.p = Popen(
            'zcat %s' % path,
            bufsize=1024,
            close_fds=True,
            shell=True,
            stdout=PIPE,
            stderr=None,
            )
        self.status = 'open'

    def __iter__(self):
        assert self.status == 'open'
        for i in self.p.stdout:
            yield i

    def next(self):
        assert self.status == 'open'
        return self.p.stdout.next()

    def close(self):
        try:
            self.p.terminate()
        except:
            from signal import SIGTERM
            try:
                os.kill(self.p.pid, SIGTERM)
            except:
                pass
        self.p.stdout.close()
        self.status = 'closed'


def sdf_iter(sdf, skip_annotations=False):

    def _next():
        if sdf == '-':
            import sys
            f = sys.stdin
        elif hasattr(sdf, 'next'):
            f = sdf
        else:
            f = file(sdf, 'r')
        buf = []
        structure_ends = False
        for i in f:
            if skip_annotations:
                if not structure_ends:
                    if i.startswith('>'):
                        structure_ends = True
                    else:
                        buf += i
                elif i.startswith('$$$$'):
                    structure_ends = False
                    buf += i
            else:
                buf += i
            if i.startswith('$$$$'):
                yield ''.join(buf)
                buf = []

    return _next()


def get_sdf_tags(sdf):
    """parse the sdf tags"""

    tag_pattern = re.compile(">\s+<([^>]+)>([^>$]+)", re.DOTALL)
    tags = tag_pattern.findall(sdf)
    tagdict = dict()

    # process each tag

    for (name, value) in tags:
        tagdict[name.strip()] = value.strip()
    return tagdict



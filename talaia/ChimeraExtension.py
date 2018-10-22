#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import chimera.extension
from Midas.midas_text import addCommand, doExtensionFunc
from talaia import enable, disable


def cmd_enable(cmdName, args):
    doExtensionFunc(enable, args, specInfo=[("spec", "selection", "residues")])


def cmd_disable(cmdName, args):
    doExtensionFunc(disable, args)


addCommand("talaia", cmd_enable, revFunc=cmd_disable)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
TALAIA
======

Simplistic 3D dictionary of amino acids using geometric shapes
for UCSF Chimera.
"""
import argparse
from time import sleep as sleep
from collections import defaultdict
from textwrap import dedent
import math
import numpy as np
import chimera
from Bld2VRML import openFileObject as openBildFileObject
from chimera import runCommand as run, cross, Point, Vector, preferences
from chimera import selection, specifier, Xform
from Shape import shapecmd
import Matrix as M
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO


AMINOACIDS = {
    "GLY": dict(
        name="GLY",
        figure="sphere",
        size=0.9,
        color1=(0, 1, 0),
        all_atoms=["N", "CA", "O", "C"],
        center="CA",
    ),
    "ALA": dict(
        name="ALA",
        figure="sphere",
        size=1.3,
        color1=(0, 1, 0),
        all_atoms=["N", "CA", "O", "C", "CB"],
        center="CA",
    ),
    "VAL": dict(
        name="VAL",
        figure="ellipsoid",
        color1=(0, 1, 0),
        all_atoms=["N", "CA", "C", "O", "CB", "CG1", "CG2"],
        k="CA",
        n="CB",
    ),
    "LEU": dict(
        name="LEU",
        figure="ellipsoid",
        color1=(0, 1, 0),
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
        k="CA",
        n="CG",
    ),
    "ILE": dict(
        name="ILE",
        figure="ellipsoid",
        color1=(0, 1, 0),
        all_atoms=["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
        k="CA",
        n="CG1",
    ),
    "MET": dict(
        name="MET",
        figure="ellipsoid",
        color1=(0, 1, 0),
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
        k="CB",
        n="SD",
    ),
    "PHE": dict(
        name="PHE",
        figure="hexagon",
        size=3.4,
        color1="green",
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        center="CG",
        k="CD2",
        n="CD1",
        p="CE2",
        u="CZ",
    ),
    "TRP": dict(
        name="TRP",
        figure="hexagon",
        size=3.5,
        color1="green",
        all_atoms=[
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "CG",
            "CD1",
            "CD2",
            "NE1",
            "CE2",
            "CE3",
            "CZ2",
            "CZ3",
            "CH2",
        ],
        center="CD1",
        k="CE2",
        n="CD2",
        p="CH2",
        u="CZ3",
        rot_factor=1.5,
    ),
    "SER": dict(
        name="SER",
        figure="cone",
        size=3.4,
        color1="aquamarine",
        color2="dark red",
        all_atoms=["N", "CA", "O", "C", "CB", "OG"],
        center="OG",
        n="CA",
        p="O",
    ),
    "THR": dict(
        name="THR",
        figure="cone",
        size=3.4,
        color1="aquamarine",
        color2="dark red",
        all_atoms=["N", "CA", "O", "C", "CB", "OG1", "CG2"],
        center="OG1",
        n="CA",
        p="O",
    ),
    "CYS": dict(
        name="CYS",
        figure="cone",
        size=3.4,
        color1="aquamarine",
        color2="yellow",
        all_atoms=["N", "CA", "C", "O", "CB", "SG"],
        center="SG",
        n="CA",
        p="O",
    ),
    "ASN": dict(
        name="ASN",
        figure="triangle",
        size=3.4,
        color1="aquamarine",
        color2=("blue", "dark red"),
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
        k="CA",
        n="CG",
        p="OD1",
        u="ND2",
    ),
    "GLN": dict(
        name="GLN",
        figure="triangle",
        size=3.4,
        color1="aquamarine",
        color2=("blue", "dark red"),
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
        k="CB",
        n="CD",
        p="OE1",
        u="NE2",
    ),
    "TYR": dict(
        name="TYR",
        figure="hexagon",
        size=3.4,
        color1=("aquamarine", "dark red"),
        all_atoms=[
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "CG",
            "CD1",
            "CD2",
            "CE1",
            "CE2",
            "CZ",
            "OH",
        ],
        center="CG",
        k="CD2",
        n="CD1",
        p="CE2",
        u="CZ",
    ),
    "LYS": dict(
        name="LYS",
        figure="rectangle",
        size=2.1,
        color1="dodger blue",
        color2="blue",
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
        center="CD",
        n="CA",
        p="O",
        k="NZ",
    ),
    "ARG": dict(
        name="ARG",
        figure="rectangle",
        size=2.2,
        color1="dodger blue",
        color2="blue",
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
        center="NE",
        n="NH1",
        p="NH2",
        k="CZ",
    ),
    "ASP": dict(
        name="ASP",
        figure="triangle",
        size=3.4,
        color1="red",
        color2=("dark red", "dark red"),
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
        k="CA",
        n="CG",
        p="OD1",
        u="OD2",
    ),
    "GLU": dict(
        name="GLU",
        figure="triangle",
        size=3.4,
        color1="red",
        color2=("dark red", "dark red"),
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
        k="CB",
        n="CD",
        p="OE1",
        u="OE2",
    ),
    "PRO": dict(
        name="PRO",
        figure="cube",
        size=1.7,
        color1="white",
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD"],
        center="N",
        n="CA",
        p="O",
    ),
    "HIS": dict(
        name="HIS",
        figure="pentagon",
        size=3.4,
        color1=["blue", "aquamarine", "dark red"],
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
        center="CG",
        n="CB",
        k="CD2",
        p="ND1",
        rot_factor=-1,
    ),
}

ALTERNATIVE_NAME = {"CYX": "CYS", "HIE": "HIS", "HE1": "HIS", "CS1": "CYS", "HID": "HIS", "HIC": "HIS"}

IONS_AND_METALS = ["MG", "NI", "PT1", "K", "CA", "NA", "F"]

all_refreshing = []


class BalloonPopup:
    """
    CLass to add popup custom balloon for vrml models.
    """
    def __init__(self, info):
        """
        info: value to show on popup balloon.
        """
        self.info = info
    
    def oslIdent(self):
        """
        return popup balloon message
        """
        return self.info



def extract_residue_message(res):
    """
    Create custom message for popup balloon for vrml models.
    """
    chain = res.id.chainId
    res_type = res.type
    pos = res.id.position
    message = res.type + " " + str(pos) + "." + chain
    return message


def get_star_atom_names(res):
    """
    get center, n and p atom names to build star for non-standard residues
    """
    clist = sorted([atom for atom in res.atoms if atom.name.startswith('C') and atom.name != 'C'])
    c = clist[0]
    n = c.neighbors[0]
    p = c.neighbors[1]
    return c, n, p


def enable(selection=None, transparency=0):
    """
    Depiction of every residue contained in selection with a vrml model.
    Parameters used for figures' depiction are fetched from AMINOACIDS dictionary.
    Parameters
    ----------
    selection: str
    	Chimera specification query or list of residues
    transparency: float
    """
    if selection is None:
        nearRes = (
            chimera.specifier.evalSpec("ligand zr < 8").residues()
            or chimera.specifier.evalSpec("all").residues())

    elif isinstance(selection, basestring):
        nearRes = chimera.specifier.evalSpec(selection).residues()
    else:  # we assume it is a list of residues. If list of atoms, residues are depicted anyway.
        nearRes = selection
    complete_residues(nearRes)
    for res in nearRes:
        message = extract_residue_message(res)
        res_info = AMINOACIDS.get(res.type)
        if res_info is None:
            if res.type in ALTERNATIVE_NAME.keys():
                alt_name = ALTERNATIVE_NAME.get(res.type)
                res_info = AMINOACIDS.get(alt_name)
            elif res.isMetal:
                continue
            elif res.type not in AMINOACIDS.keys() and not res.isIsolated:
                name = res.type
                center, n, p = get_star_atom_names(res)
                size = 3
                color = "orange"  # "dim gray"
                print(center.coord(), n.coord(), p.coord(), size, color, name, message, transparency)
                star(center.coord(), n.coord(), p.coord(), size, color, name, message, transparency)         
            else:
                if "CA" not in res.atomsMap:
                    center = res.atoms[0].coord()
                    continue
        if res.type == 'HIS':
            rt = proton_eval(res)
        rt = proton_eval(res)

        if res_info:
            res._depicted = True
            res.molecule._depicted = dict(
            	   selection=selection, transparency=transparency
            )
            name = res_info.get("name")
            shape = res_info.get("figure")
            size = res_info.get("size")
            r1 = size
            tp = transparency
            color1 = res_info.get("color1")
            color2 = res_info.get("color2")

            if res_info.get("center") is not None:
                center = res.atomsMap[res_info["center"]][0].coord()
                q = preparation(center)
            if res_info.get("n") is not None:
                n = res.atomsMap[res_info["n"]][0].coord()
            if res_info.get("p") is not None:
                p = res.atomsMap[res_info["p"]][0].coord()
            if res_info.get("k") is not None:
                k = res.atomsMap[res_info["k"]][0].coord()
            if res_info.get("u") is not None:
                u = res.atomsMap[res_info["u"]][0].coord()

            if shape == "sphere":
                sph = shapecmd.sphere_shape(radius=r1, center=q, color=(0, 1, 0, 1 - tp))
                sph.molecule = BalloonPopup("#? " + message)
            elif shape == "ellipsoid":
                ellipsoid(k, n, message, tp)
            elif shape == "triangle":
                triangle(n, p, k, u, r1, color1, color2, tp, name, message)
            elif shape == "cone":
                cone(center, n, p, r1, color1, color2, name, message, tp)
            elif shape == "cube":
                cube(center, n, p, r1, color1, name, message, tp)
            elif shape == "rectangle":
                rectangle(center, n, p, k, r1, color1, color2, name, message, tp)
            elif shape == "pentagon":
                pentagon(center, n, p, k, r1, color1, rt, name, message, tp)
            elif shape == "hexagon":
                hexagon(center, n, p, k, u, r1, color1, name, message, tp)
            
            opened = chimera.openModels.list()
            key = 'refreshing_' + str(res.id)
            key = dict(
                residue=res,
                name=name,
                shape=shape,
                frames=dict(ref_center=detect_type(res, name)[0],
                            ref_points=detect_type(res, name)[1],
                            vec2=None)
            )
            # As MET has 2 models associated, both need to be added to all_refreshing so MD works properly
            if name == "MET":
                key_met = key.copy()
                key_met["model"]=opened[-2]
                all_refreshing.append(key_met)
            key["model"] = opened[-1]
            all_refreshing.append(key)

    subscribe_events()
    return all_refreshing


def disable():
    """
    Removes all vrml models on display.
    """
    all_models = chimera.openModels.list()
    for item in all_models:
        if (
                item.name.startswith("vrml_")
                or item.name == "ellipsoid"
                or item.name == "sphere"
                or item.name == "cylinder"
        ):
            chimera.openModels.close([item])
        if hasattr(item, "_depicted"):
            del item._depicted
            for r in item.residues:
                if hasattr(r, "_depicted"):
                    del r._depicted

    unsubscribe_events()
    for d in all_refreshing:
        d.clear()
    all_refreshing[:] = []


def subscribe_events():
    """
    """
    chimera._depicter_handler = chimera.triggers.addHandler("Molecule", callback, all_refreshing)


def unsubscribe_events():
    """
    """
    if hasattr(chimera, "_depicter_handler"):
        chimera.triggers.deleteHandler("Molecule", chimera._depicter_handler)
        del chimera._depicter_handler


def callback(name, data, changes):
    """
    Update shapes position and orientation after coordinates change.
    """
    depicted = [m for m in chimera.openModels.list() if hasattr(m, "_depicted")]
    options = depicted[0]._depicted if depicted else {}
    modified_mols = set(depicted) & changes.modified
    deleted_mols = set(depicted) & changes.deleted
    if any(
            [
                (modified_mols and "activeCoordSet changed" in changes.reasons),
                (modified_mols and "atoms moved" in changes.reasons),
                deleted_mols]
    ):
        for i in range(len(all_refreshing)):
            ref_points = all_refreshing[i]['frames'].get('ref_points')
            ref_center = all_refreshing[i]['frames'].get('ref_center')
            return_coord, s2_points = detect_type(all_refreshing[i].get('residue'),
                                                  all_refreshing[i].get('name'))
            refresh(all_refreshing[i]['shape'], ref_center, return_coord, ref_points,
                    s2_points, all_refreshing[i].get('model'))


def proton_eval(residue):
    """
    Identify protonation state of Histidine residues.
    Returns value between 0 and 2, used to determine colour
    of histidine pentagon according to protonation.
    """
    ct = 0
    h_ct = 0
    for key in residue.atomsMap.keys():
        k = str(key)
        if k.startswith('H'):
            h_ct += 1
    if h_ct != 0:
        for at in residue.atoms:
            if at.idatmType == 'N2':
                ct += 1
    else:
        ct = 1

    return ct


def Zmatrix(origin_coord):
    """
    """
    ZERO = Point(0, 0, 0)
    get_zero = ZERO - origin_coord
    zero_matrix = Xform.translation(get_zero)
    return zero_matrix


def Tmatrix(return_coord):
    """
    """
    ZERO = Point(0, 0, 0)
    get_back = return_coord - ZERO
    translation_matrix = Xform.translation(get_back)
    return translation_matrix


def Rmatrix(vec1, vec2):
    """
    """
    vec1.normalize()
    vec2.normalize()
    near_zero = 1e-15
    rvec = chimera.cross(vec1, vec2)
    if rvec == chimera.Vector(0, 0, 0):
        MI = ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0))
        return M.chimera_xform(MI)
    else:
        theta = chimera.angle(vec1, vec2)
        rotation_matrix = Xform.rotation(rvec, theta)
        return rotation_matrix


def refresh(shape, origin_coord, return_coord, ref_points, s2_points, model):
    """
    """
    if shape == 'sphere':
        xform = sphere_transform(origin_coord, return_coord)
    else:
        xform = transform(origin_coord, return_coord, ref_points, s2_points)
    mol = chimera.openModels.list()[0]
    curxform = mol.openState.xform
    XF = M.xform_matrix(xform)
    CX = M.xform_matrix(curxform)
    mlist = [CX, XF]
    to_apply = M.chimera_xform(M.multiply_matrices(*mlist))
    model.openState.xform = to_apply


def key_coords(res, a1, a2=None, a3=None):
    """
    """
    if a2:
        if a3:
            for at in res.atoms:
                if at.name == a1:
                    atom1 = at.coord()
                elif at.name == a2:
                    atom2 = at.coord()
                elif at.name == a3:
                    atom3 = at.coord()
            return atom1, atom2, atom3
        else:
            for at in res.atoms:
                if at.name == a1:
                    atom1 = at.coord()
                elif at.name == a2:
                    atom2 = at.coord()
            return atom1, atom2
    else:
        for at in res.atoms:
            if at.name == a1:
                atom1 = at.coord()
        return atom1


def detect_non_standard(res):
    """
    """
    atom_map = res.atomsMap.keys()
    for item in AMINOACIDS:
        calls = {}
        if compare_all_atoms(AMINOACIDS[item], atom_map, calls):
            return res, item


def compare_all_atoms(pattern, text, calls):
    """
    """
    text = convert_unicode_to_string(remove_hydrogens(text))
    key = str(pattern) + ':' + str(text)
    if key in calls:
        return calls[key]
    if len(pattern) == 0:
        return len(text)
    elif len(text) == 0:
        return len(pattern)
    else:
        if pattern[0] == text[0]:
            m_cost = compare_all_atoms(pattern[1:], text[1:], calls)
        else:
            m_cost = compare_all_atoms(pattern[1:], text[1:], calls) + 1
        i_cost = compare_all_atoms(pattern[:], text[1:], calls) + 1
        d_cost = compare_all_atoms(pattern[1:], text[:], calls)
        minimum = min(m_cost, i_cost, d_cost)
        calls[str(pattern) + ':' + str(text)] = minimum
        return minimum


def convert_unicode_to_string(map_list):
    """
    """
    new_list = []
    for i in map_list:
        j = str(i)
        new_list.append(j)
    return new_list


def remove_hydrogens(atom_list):
    """
    """
    no_hydrogen = []
    for i in atom_list:
        if i.startswith('H'):
            pass
        else:
            no_hydrogen.append(i)
    return no_hydrogen


def detect_type(res, res_name):
    """
    non_standard = []
    if res.type not in AMINOACIDS.keys():
    	non_standard.append(res)

    """
    if res_name == 'GLY' or res_name == 'ALA':
        center = key_coords(res, AMINOACIDS['GLY']['center'])
        return center, [center]
    elif res_name == 'VAL':
        k, n = key_coords(res, AMINOACIDS['VAL']['k'], AMINOACIDS['VAL']['n'])
        center = k + ((n - k) / 2)
        return center, [center, n, k]
    elif res_name == 'LEU':
        k, n = key_coords(res, AMINOACIDS['LEU']['k'], AMINOACIDS['LEU']['n'])
        center = k + ((n - k) / 2)
        return center, [center, n, k]
    elif res_name == 'ILE':
        k, n = key_coords(res, AMINOACIDS['ILE']['k'], AMINOACIDS['ILE']['n'])
        center = k + ((n - k) / 2)
        return center, [center, n, k]
    elif res_name == 'MET':
        k, n = key_coords(res, AMINOACIDS['MET']['k'], AMINOACIDS['MET']['n'])
        center = k + ((n - k) / 2)
        return center, [center, center, n]
    elif res_name == 'SER':
        center, n = key_coords(res, AMINOACIDS['SER']['n'], AMINOACIDS['SER']['center'])
        return center, [center, center, n]
    elif res_name == 'THR':
        center, n = key_coords(res, AMINOACIDS['THR']['n'], AMINOACIDS['THR']['center'])
        return center, [center, center, n]
    elif res_name == 'CYS':
        center, n = key_coords(res, AMINOACIDS['CYS']['n'], AMINOACIDS['CYS']['center'])
        return center, [center, center, n]
    elif res_name == 'ASP':
        n, k = key_coords(res, AMINOACIDS['ASP']['n'], AMINOACIDS['ASP']['k'])
        p, u = key_coords(res, AMINOACIDS['ASP']['p'], AMINOACIDS['ASP']['u'])
        s2_points = calculate_triangle_vertex(n, k, p, u)
        return s2_points[0], s2_points
    elif res_name == 'ASN':
        n, k = key_coords(res, AMINOACIDS['ASN']['n'], AMINOACIDS['ASN']['k'])
        p, u = key_coords(res, AMINOACIDS['ASN']['p'], AMINOACIDS['ASN']['u'])
        s2_points = calculate_triangle_vertex(n, k, p, u)
        return s2_points[0], s2_points
    elif res_name == 'GLU':
        n, k = key_coords(res, AMINOACIDS['GLU']['n'], AMINOACIDS['GLU']['k'])
        p, u = key_coords(res, AMINOACIDS['GLU']['p'], AMINOACIDS['GLU']['u'])
        s2_points = calculate_triangle_vertex(n, k, p, u)
        return s2_points[0], s2_points
    elif res_name == 'GLN':
        n, k = key_coords(res, AMINOACIDS['GLN']['n'], AMINOACIDS['GLN']['k'])
        p, u = key_coords(res, AMINOACIDS['GLN']['p'], AMINOACIDS['GLN']['u'])
        s2_points = calculate_triangle_vertex(n, k, p, u)
        return s2_points[0], s2_points
    elif res_name == 'ARG':
        center, n = key_coords(res, AMINOACIDS['ARG']['center'], AMINOACIDS['ARG']['n'])
        p, k = key_coords(res, AMINOACIDS['ARG']['p'], AMINOACIDS['ARG']['k'])
        return center, [center, n, p, center, k]
    elif res_name == 'LYS':
        center, n, k = key_coords(res, AMINOACIDS['LYS']['center'], AMINOACIDS['LYS']['n'],
                                  AMINOACIDS['LYS']['k'])
        return center, [center, n, k]
    elif res_name == 'TRP':
        center, k, n = key_coords(res, AMINOACIDS['TRP']['center'], AMINOACIDS['TRP']['k'],
                                  AMINOACIDS['TRP']['n'])
        p, u = key_coords(res, AMINOACIDS['TRP']['p'], AMINOACIDS['TRP']['u'])
        fcenter = ((k - p) * (center.distance(u) / ((-2)*k.distance(p)))) + center
        return fcenter, [center, k, n, k, p]
    elif res_name == 'TYR':
        center, k, n = key_coords(res, AMINOACIDS['TYR']['center'], AMINOACIDS['TYR']['k'],
                                  AMINOACIDS['TYR']['n'])
        p, u = key_coords(res, AMINOACIDS['TYR']['p'], AMINOACIDS['TYR']['u'])
        fcenter = ((k - p) * (center.distance(u) / ((-2)*k.distance(p)))) + center
        return fcenter, [center, k, n, k, p]
    elif res_name == 'PHE':
        center, k, n = key_coords(res, AMINOACIDS['TYR']['center'], AMINOACIDS['TYR']['k'],
                                  AMINOACIDS['TYR']['n'])
        p, u = key_coords(res, AMINOACIDS['TYR']['p'], AMINOACIDS['TYR']['u'])
        fcenter = ((k - p) * (center.distance(u) / ((-2)*k.distance(p)))) + center
        return fcenter, [center, k, n, k, p]
    elif res_name == 'HIS':
        center, n = key_coords(res, AMINOACIDS['HIS']['center'], AMINOACIDS['HIS']['n'])
        k, p = key_coords(res, AMINOACIDS['HIS']['k'], AMINOACIDS['HIS']['p'])
        center2 = ((-1)*(n - center)) + center
        return center2, [center, n, center2, k, p]
    elif res_name == 'PRO':
        center, n = key_coords(res, AMINOACIDS['PRO']['center'], AMINOACIDS['PRO']['n'])
        return center, [center, n, center]
    else:
        print 'residue no detectat'


def sphere_transform(origin_coord, return_coord):
    """
    """
    zero_matrix = Zmatrix(origin_coord)
    MZ = M.xform_matrix(zero_matrix)

    translation_matrix = Tmatrix(return_coord)
    MT = M.xform_matrix(translation_matrix)

    mlist = [MT, MZ]
    to_apply = M.chimera_xform(M.multiply_matrices(*mlist))

    return to_apply


def transform(origin_coord, return_coord, ref_points, s2_points):
    """
    """
    zero_matrix = Zmatrix(origin_coord)
    MZ = M.xform_matrix(zero_matrix)

    translation_matrix = Tmatrix(return_coord)
    MT = M.xform_matrix(translation_matrix)

    vec1 = ref_points[1] - ref_points[2]
    vec2 = s2_points[1] - s2_points[2]
    rotation_matrix = Rmatrix(vec1, vec2)
    MR = M.xform_matrix(rotation_matrix)

    mlist = [MT, MR, MZ]
    MA = M.multiply_matrices(*mlist)
    to_apply = M.chimera_xform(MA)

    if len(ref_points) == 5:
        x_points = recalculate_vertex(to_apply, ref_points)
        vecx = x_points[3] - x_points[4]
        vec3 = s2_points[3] - s2_points[4]
        rotation_matrix2 = Rmatrix(vecx, vec3)
        MR2 = M.xform_matrix(rotation_matrix2)

        mlist2 = [MT, MR2, MR, MZ]
        to_apply2 = M.chimera_xform(M.multiply_matrices(*mlist2))
        return to_apply2
    else:
        return to_apply


def calculate_triangle_vertex(n, k, p, u):
    """
    """
    vec_nk = n - k
    vec_pu = p - u
    half_length = 3.4 / 2.0
    thickness = 3.4 / 6.0
    adjustment_nk = half_length / k.distance(n)
    adjustment_pu = half_length / p.distance(u)
    adj_vec_nk1 = vec_nk * adjustment_nk
    adj_vec_nk2 = adj_vec_nk1 * -1
    center = ((adj_vec_nk2) * 0.5) + n
    adj_vec_pu1 = vec_pu * (adjustment_pu)
    adj_vec_pu2 = adj_vec_pu1 * -1

    x1 = adj_vec_nk1 + center
    x2 = adj_vec_nk2 + center
    o1 = adj_vec_pu1 + x1
    o2 = adj_vec_pu1 + x2
    o3 = adj_vec_pu2 + x1

    perp_for = normalize(cross(o1 - o2, o3 - o1)) * thickness
    perp_back = -1 * perp_for

    s1 = perp_for + x2
    s2 = perp_back + x2
    s3 = perp_for + o3
    s4 = perp_back + o3
    s5 = perp_for + o1
    s6 = perp_back + o1

    return [center, s3, s5, s4, s3]


def recalculate_vertex(matrix, s1_points):
    """
    Recalculate vertex of triangle after rotation.
    """
    new_points = []
    for i in s1_points:
        new_points.append(matrix.apply(i))
    return new_points


def _vmd_trans_angle(a, b, c, delta):
    """
    Simulates VMD's `trans angle` command
    """
    ZERO = Point(0, 0, 0)
    xf = chimera.Xform.translation(b - ZERO)
    xf.rotate(cross(a - b, b - c), delta)
    xf.translate(ZERO - b)
    return xf


def _rotate(a, b, c, delta, x):
    """
    """
    rotation_xform = _vmd_trans_angle(a, b, c, delta)
    rotation_array = np.array(rotation_xform.getOpenGLMatrix()).reshape(4, 4).T
    return Vector(*np.dot(rotation_array, x.data() + (0,))[:3])


def normalize(v):
    return Vector(*M.normalize_vector(v))


def preparation(center):
    """
    Correct coordinates string to be used by shapecmd
    """
    return str(center).replace(" ", ",")


def vectorial_product(v1, v2):
    """
    Obtain product vector from two other vectors.
    """
    x = (v1.y * v2.z) - (v1.z * v2.y)
    y = (v1.x * v2.z) - (v1.z * v2.x)
    z = (v1.x * v2.y) - (v1.y * v2.x)
    return (x, y, z)


def length(v):
    """
    Obtain length of a vector
    """
    return math.sqrt(dotproduct(v, v))


def dotproduct(v1, v2):
    """
    Obtain dot product of two vectors
    """
    return sum((a * b) for a, b in zip(v1, v2))


def complete_residues(residues):
    """
    Detect any residue with missing atoms and substitutes it
    with the best rotamer from Dunbrack 2010 library.
    """
    import Rotamers

    for res in residues:
        if res.isHet:
            continue
        res_incomplete = []
        for item in AMINOACIDS:
            if res.type == AMINOACIDS[item]["name"]:
                info = AMINOACIDS[item]
                atoms = info.get("all_atoms")
                if len(res.atoms) < len(atoms):
                    res_incomplete.append(res)
                Rotamers.useBestRotamers(res.type, res_incomplete)


def ellipsoid(k, n, message, transparency=0):
    """
    Used for residues V, L, I, M.
    Define ellipsoids orientation following residue's side chains.
    """
    # calculate center coordinates with vector between k and n points.
    k = Point(*k)
    n = Point(*n)
    vec_kn = k - n
    center = (vec_kn * 0.5) + n
    q = preparation(center)
    tp = 1 - transparency
    xvec = Vector(1, 0, 0)

    x, y, z, theta = shapecmd_rotation(vec_kn, xvec)
    
    ellip = shapecmd.sphere_shape(
            radius=(2, 1, 1), center=q, rotation=(x, y, z, theta),
            color=(0, 1, 0, tp)
    )
    ellip.molecule = BalloonPopup("#? " + message)
    if "MET" in message:
        zvec = Vector(0, 0, 1)
        x2, y2, z2, theta2 = shapecmd_rotation(vec_kn, zvec)
        center = preparation(n)
        cyl = shapecmd.cylinder_shape(radius=0.9, height=0.7, slab=0.1,
                                      center=center,
                                      rotation=(x2,y2,z2,theta2),
                                      color=(1,1,0,tp))
        cyl.molecule = BalloonPopup("#? " + message)      



def ellipsoid_old(k, n, message, transparency=0):
    """
    Used for residues V, L, I, M.
    Define ellipsoids orientation following residue's side chains.
    """
    # calculate center coordinates with vector between k and n points.
    k = Point(*k)
    n = Point(*n)
    vec_kn = k - n
    center = (vec_kn * 0.5) + n
    q = preparation(center)

    tp = 1 - transparency
    xvec = Vector(1, 0, 0)

    # calculate theta, cosinus and sinus of theta
    cos_theta = dotproduct(xvec, vec_kn) / length(vec_kn) * (length(xvec))
    sin_theta = math.sqrt(1 - (cos_theta ** 2))
    theta = math.degrees(math.acos(cos_theta))
    # calculate rotation vector and convert it to a unit vector
    x, y, z = vectorial_product(xvec, vec_kn)
    rot_vec = Vector(x, y, z)
    u_rot_vec = rot_vec * (1 / length(rot_vec))

    # Apply Rodrigues rotation formula:
    # vrot = a + b + c =
    #      = v * cos_theta + (k x v) * sin_theta + k * (k dot v) * (1 - cos_theta)
    a = xvec * cos_theta
    vp = vectorial_product(u_rot_vec, xvec)
    vproduct = Vector(vp[0], vp[1], vp[2])
    b = vproduct * sin_theta
    c = u_rot_vec * dotproduct(u_rot_vec, xvec) * (1 - cos_theta)
    vrot = a + b + c

    # Check if vrot is paralel to vec_kn: Dot product = 1
    # If dotproduct != 1 neet to change y component sign.
    u_vec_kn = vec_kn * (1 / length(vec_kn))
    if (
            abs(dotproduct(vrot, u_vec_kn)) < 0.999999
            or abs(dotproduct(vrot, u_vec_kn)) > 1.0001
    ):
        y = -1 * y
    ellip = shapecmd.sphere_shape(
            radius=(2, 1, 1), center=q, rotation=(x, y, z, theta),
            color=(0, 1, 0, tp)
    )
    ellip.molecule = BalloonPopup("#? " + message)


def shapecmd_rotation(vec, ref_vec):
    """
    purpose: to calculate rotation vector and angle for ellipsoids and cylinders
        acording to the main vector of the residue in question and the 
        reference vector necessary (ellipsoids:xvec, cylinders:zvec)
    both arguments must be chimera Vector objects.
    output is values x, y, z and theta to pass as the rotation argument for the shapecmd functions
    """
    cos_theta = dotproduct(ref_vec, vec) / length(vec) * (length(ref_vec))
    sin_theta = math.sqrt(1 - (cos_theta ** 2))
    theta = math.degrees(math.acos(cos_theta))
    # calculate rotation vector and convert it to a unit vector
    x, y, z = vectorial_product(ref_vec, vec)
    rot_vec = Vector(x, y, z)
    u_rot_vec = rot_vec * (1 / length(rot_vec))

    # Apply Rodrigues rotation formula:
    # vrot = a + b + c =
    #      = v * cos_theta + (k x v) * sin_theta + k * (k dot v) * (1 - cos_theta)
    a = ref_vec * cos_theta
    vp = vectorial_product(u_rot_vec, ref_vec)
    vproduct = Vector(vp[0], vp[1], vp[2])
    b = vproduct * sin_theta
    c = u_rot_vec * dotproduct(u_rot_vec, ref_vec) * (1 - cos_theta)
    vrot = a + b + c

    # Check if vrot is paralel to vec: Dot product = 1
    # If dotproduct != 1 neet to change y component sign.
    u_vec = vec * (1 / length(vec))
    if (
            abs(dotproduct(vrot, u_vec)) < 0.999999
            or abs(dotproduct(vrot, u_vec)) > 1.0001
    ):
        y = -1 * y
    return x, y, z, theta


def cube(center, n, p, size, color1, name, message, transparency=0):
    """
    Used for resiude P.
    Define a white cubic shape from three atom coordinates.
    """
    center = Point(*center)
    n = Point(*n)
    p = Point(*p)
    tp = transparency
    # create vector between n and center and adjust its length.
    vec_AB = n - center
    half_length = size / 2.0
    adjustment = half_length / center.distance(n)
    adj_vec_AB_1 = vec_AB * adjustment
    adj_vec_AB_2 = -1 * adj_vec_AB_1
    # obtain two points from center with vector's direction.
    # calculate perpendicual vectors from them and point p.
    x1 = adj_vec_AB_1 + center
    x2 = adj_vec_AB_2 + center
    perp_1 = normalize(cross(x1 - x2, x2 - p))
    perp1 = perp_1 * half_length
    perp2 = -1 * perp1
    # obtain 4 points describing a square on the plane defined by
    # vec_AB and perp_1
    o1 = perp1 + x1
    o2 = perp1 + x2
    o3 = perp2 + x1
    o4 = perp2 + x2
    # calculate perpendicular vector and obtain final 8 points.
    perp_for = normalize(cross(o1 - o2, o3 - o1)) * half_length
    perp_back = -1 * perp_for
    points = dict(
        s1=perp_for + o1,
        s2=perp_for + o2,
        s3=perp_for + o3,
        s4=perp_for + o4,
        s5=perp_back + o1,
        s6=perp_back + o2,
        s7=perp_back + o3,
        s8=perp_back + o4,
    )
    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {s2} {s3} {s4}
    .polygon {s1} {s2} {s6}
    .polygon {s4} {s7} {s8}
    .polygon {s5} {s6} {s8}
    .polygon {s2} {s4} {s8}
    .polygon {s1} {s5} {s7}
    .polygon {s1} {s3} {s2}
    .polygon {s3} {s7} {s4}
    .polygon {s1} {s6} {s5}
    .polygon {s5} {s8} {s7}
    .polygon {s2} {s8} {s6}
    .polygon {s1} {s7} {s3}
    """.format(
        color1=color1, tp=tp, **points
    )
    add_vrml_model(bild, "talaia " + message, message)


def triangle(n, p, k, u, size, color1, color2, transparency, name, message):
    """
    Used for residues E, D, Q, N.
    Define oriented orange triangle with one side of different
    color depending on residue type.
    """
    n = Point(*n)
    p = Point(*p)  # O
    k = Point(*k)
    u = Point(*u)  # N
    tp = transparency
    # create two vectors and adapt their length to calculate center point.
    vec_nk = n - k
    vec_pu = p - u  # amide vector
    half_length = size / 2.0
    thickness = size / 6.0
    adjustment_nk = half_length / k.distance(n)
    adjustment_pu = half_length / p.distance(u)
    adj_vec_nk1 = vec_nk * adjustment_nk
    adj_vec_nk2 = adj_vec_nk1 * -1
    center = ((adj_vec_nk2) * 0.5) + n
    adj_vec_pu1 = vec_pu * (adjustment_pu)
    adj_vec_pu2 = adj_vec_pu1 * -1

    x1 = adj_vec_nk1 + center
    x2 = adj_vec_nk2 + center
    o1 = adj_vec_pu1 + x1
    o2 = adj_vec_pu1 + x2
    o3 = adj_vec_pu2 + x1

    perp_for = normalize(cross(o1 - o2, o3 - o1)) * thickness
    perp_back = -1 * perp_for
    points = dict(
        s1=perp_for + x2,
        s2=perp_back + x2,
        s3=perp_for + o3,
        s4=perp_back + o3,
        s5=perp_for + o1,
        s6=perp_back + o1,
        d1=perp_for + x1,
        d2=perp_back + x1,
    )

    # unpack base colors
    color21 = color2[0]
    color22 = color2[1]

    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {s1} {s2} {s3}
    .polygon {s2} {s3} {s4}
    .polygon {s1} {s2} {s5}
    .polygon {s2} {s6} {s5}
    .color {color22}
    .polygon {s6} {s5} {d1}
    .polygon {s6} {d1} {d2}
    .color {color21}
    .polygon {s3} {d2} {d1}
    .polygon {s3} {s4} {d2}
    .color {color1}
    .polygon {s1} {s3} {s5}
    .polygon {s2} {s4} {s6}
    """.format(
        color1=color1, color22=color22, color21=color21, tp=tp, **points
    )
    add_vrml_model(bild, "talaia " + message, message)
    return(points)


def cone(center, n, p, size, color1, color2, name, message, transparency=0):
    """
    Used for residues S, T, C.
    Defines an orange 3-base pyramid with vertex pointing to OH/SH
    resiudes' functional group to show its orientation.
    """
    n = Point(*n)
    center = Point(*center)
    p = Point(*p)
    tp = transparency

    vec_AB = n - center
    half_length = size / 2.0

    _adjustment = half_length / center.distance(n)
    adjustment1 = _adjustment * 0.66
    adjustment2 = _adjustment * 1.33

    adj_vec_AB_1 = adjustment1 * vec_AB
    adj_vec_AB_2 = -adjustment2 * vec_AB

    x1 = adj_vec_AB_1 + center
    x2 = adj_vec_AB_2 + center

    perp1 = normalize(cross(x1 - x2, x2 - p))
    o1 = perp1 + x1
    perp3 = normalize(cross(o1 - x2, x2 - n)) * half_length
    o3 = perp3 + x1

    d1 = o3 - x1

    vertex = [_rotate(o1, o3, x1, alpha, d1) + x1 for alpha in (120, 240, 360)]

    x4 = ((x2.distance(center)) / (length(vec_AB))) * vec_AB + x2

    outer = []

    for v in vertex:
        outer.append((((x2.distance(center)) / (length(vec_AB))) * vec_AB) + v)

    color_and_points = {
        "outer_{}".format(i + 1): value for i, value in enumerate(outer)
    }

    mid_list = []
    for item in color_and_points:
        vec = (color_and_points[item] - x4) / 1.5
        p = x4 + vec
        mid_list.append(p)

    color_and_points.update(dict(mid1=mid_list[0], mid2=mid_list[1],
                                 mid3=mid_list[2], x1=x1, x4=x4, color1=color1,
                                 color2=color2, tp=tp))
    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {mid3} {outer_1} {outer_2}
    .polygon {mid3} {mid2} {outer_2}
    .polygon {mid2} {outer_2} {outer_3}
    .polygon {mid2} {mid1} {outer_3}
    .polygon {mid1} {outer_3} {outer_1}
    .polygon {mid1} {mid3} {outer_1}
    .polygon {outer_1} {outer_2} {outer_3}
    .color {color2}
    .transparency {tp}
    .polygon {x4} {mid1} {mid2}
    .polygon {x4} {mid2} {mid3}
    .polygon {x4} {mid1} {mid3}
    """.format(**color_and_points)

    add_vrml_model(bild, "talaia " + message, message)


def rectangle(center, n, p, k, size, color1, color2, name, message, transparency=0):
    """
    Used for residues K and R.
    Describe a blue rectangle placed on the end of the residue side chain
    with a side coloured with a darker blue to indicate orientation.
    """
    center = Point(*center)
    n = Point(*n)
    p = Point(*p)
    k = Point(*k)
    tp = transparency
    thickness = size / 1.7
    half_length = size / 1.2

    # Choose between K (size 2.1) and R (size 2.2)
    if size == 2.1:
        vec_AB = n - k
        adjustment = half_length / k.distance(n)
        adj_vec_AB_1 = adjustment * vec_AB
        adj_vec_AB_2 = -adjustment * vec_AB

        x1 = adj_vec_AB_1 + center
        x2 = adj_vec_AB_2 + center

        perp1 = normalize(cross(x1 - x2, x2 - p)) * half_length
        perp2 = perp1 * -1
    elif size == 2.2:
        vec_np = n - p
        vec_kc = center - k
        adjust_length = half_length / center.distance(k)
        adj_vec_kc_1 = adjust_length * vec_kc
        adj_vec_kc_2 = -1 * adj_vec_kc_1

        x1 = adj_vec_kc_1 + center
        x2 = adj_vec_kc_2 + center

        perp1 = normalize(cross(vec_np, vec_kc)) * thickness
        perp2 = -1 * perp1

    o1 = perp1 + x1
    o2 = perp1 + x2
    o3 = perp2 + x1
    d1 = x2 - center

    outer = [
        _rotate(o1, center, o2, alpha, d1) + center for alpha in (22, 157, 202, 337)
    ]

    perp_for = thickness * normalize(cross(o1 - o2, o3 - o1))
    perp_back = perp_for * -1
    color_and_points = dict(
        center_1=perp_for + center,
        center_2=perp_back + center,
        front_1=perp_for + outer[0],
        front_2=perp_for + outer[1],
        front_3=perp_for + outer[2],
        front_4=perp_for + outer[3],
        back_1=perp_back + outer[0],
        back_2=perp_back + outer[1],
        back_3=perp_back + outer[2],
        back_4=perp_back + outer[3],
        color1=color1,
        color2=color2,
        tp=tp,
    )
    # Draw the rectangle
    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {front_1} {front_2} {front_3}
    .polygon {front_1} {front_3} {front_4}
    .polygon {back_1} {back_3} {back_2}
    .polygon {back_1} {back_4} {back_3}
    .polygon {back_1} {back_2} {front_1}
    .polygon {back_2} {front_2} {front_1}
    .polygon {back_2} {back_3} {front_2}
    .polygon {back_3} {front_3} {front_2}
    .polygon {back_3} {back_4} {front_3}
    .polygon {back_4} {front_4} {front_3}
    .color {color2}
    .polygon {back_4} {back_1} {front_4}
    .polygon {back_1} {front_1} {front_4}
    """.format(
        **color_and_points
    )
    add_vrml_model(bild, "talaia " + message, message)


def ion_star(center, size, color1, transparency=0):
    """
    Used for metal atoms and ions.
    Describes a grey 3D star centered on the only atom of the residue.
    """
    center = Point(*center)
    zero = Point(0, 0, 0)
    v0 = center - zero
    nv0 = v0 * (1 / length(v0))
    xref = Vector(1, 0, 0)
    n = (size * nv0) + center
    x, y, z = vectorial_product(nv0, xref)
    vt = Vector(x, y, z)
    p = (size * vt) + center
    tp = transparency

    vec_AB = n - center
    shape_size = size * 1.5
    half_length = shape_size / 2.0
    thickness = shape_size / 4.0

    adjustment = half_length / center.distance(n)
    adj_vec_AB_1 = adjustment * vec_AB
    adj_vec_AB_2 = -1 * adj_vec_AB_1

    x1 = adj_vec_AB_1 + center
    x2 = adj_vec_AB_2 + center

    perp1 = normalize(cross(x1 - x2, x2 - p)) * half_length
    perp2 = -1 * perp1

    o1 = perp1 + x1
    o2 = perp1 + x2
    o3 = perp2 + x1
    d1 = x2 - center
    d2 = (x1 - center) * 0.5

    # Rotate by 72 degrees to identify other points of the star
    outer = [
        _rotate(o1, center, o2, alpha, d1) + center
        for alpha in (72, 144, 216, 288, 360)
    ]
    inner = [
        _rotate(o1, center, o2, alpha, d2) + center
        for alpha in (72, 144, 216, 288, 360)
    ]

    perp_for = thickness * normalize(cross(o1 - o2, o3 - o1))
    perp_back = -1 * perp_for
    center_1 = perp_for + center
    center_2 = perp_back + center

    color_and_points = {
        "outer_{}".format(i + 1): value for i, value in enumerate(outer)
    }
    color_and_points.update(
        {"inner_{}".format(i + 1): value for i, value in enumerate(inner)}
    )
    color_and_points.update(
        dict(color1=color1, tp=tp, center_1=center_1, center_2=center_2)
    )

    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {outer_1} {center_1} {inner_3}
    .polygon {outer_1} {inner_3} {center_2}
    .polygon {outer_1} {inner_4} {center_1}
    .polygon {outer_1} {center_2} {inner_4}
    .polygon {outer_2} {center_1} {inner_4}
    .polygon {outer_2} {inner_4} {center_2}
    .polygon {outer_2} {inner_5} {center_1}
    .polygon {outer_2} {center_2} {inner_5}
    .polygon {outer_3} {center_1} {inner_5}
    .polygon {outer_3} {inner_5} {center_2}
    .polygon {outer_3} {inner_1} {center_1}
    .polygon {outer_3} {center_2} {inner_1}
    .polygon {outer_4} {center_1} {inner_1}
    .polygon {outer_4} {inner_1} {center_2}
    .polygon {outer_4} {inner_2} {center_1}
    .polygon {outer_4} {center_2} {inner_2}
    .polygon {outer_5} {center_1} {inner_2}
    .polygon {outer_5} {inner_2} {center_2}
    .polygon {outer_5} {inner_3} {center_1}
    .polygon {outer_5} {center_2} {inner_3}
    """.format(
        **color_and_points
    )
    add_vrml_model(bild, "ion_star")


def star(center, n, p, size, color1, name, message, transparency=0):
    """
    Used for residues out of AMINOACIDS dictionary and ALTERNATIVE_NAMES
    dictionary that are not HETATM (solvent, metal, ligand, cofactor).
    Describes a red 3D star centered on CA atom.
    """
    n = Point(*n)
    center = Point(*center)
    p = Point(*p)
    tp = transparency

    vec_AB = n - center
    shape_size = size * 1.5
    half_length = shape_size / 2.0
    thickness = shape_size / 4.0

    adjustment = half_length / center.distance(n)
    adj_vec_AB_1 = adjustment * vec_AB
    adj_vec_AB_2 = -1 * adj_vec_AB_1

    x1 = adj_vec_AB_1 + center
    x2 = adj_vec_AB_2 + center

    perp1 = normalize(cross(x1 - x2, x2 - p)) * half_length
    perp2 = -1 * perp1

    o1 = perp1 + x1
    o2 = perp1 + x2
    o3 = perp2 + x1
    d1 = x2 - center
    d2 = (x1 - center) * 0.5

    # Rotate by 72 degrees to identify other points of the star
    outer = [
        _rotate(o1, center, o2, alpha, d1) + center
        for alpha in (72, 144, 216, 288, 360)
    ]
    inner = [
        _rotate(o1, center, o2, alpha, d2) + center
        for alpha in (72, 144, 216, 288, 360)
    ]

    perp_for = thickness * normalize(cross(o1 - o2, o3 - o1))
    perp_back = -1 * perp_for
    center_1 = perp_for + center
    center_2 = perp_back + center

    color_and_points = {
        "outer_{}".format(i + 1): value for i, value in enumerate(outer)
    }
    color_and_points.update(
        {"inner_{}".format(i + 1): value for i, value in enumerate(inner)}
    )
    color_and_points.update(
        dict(color1=color1, tp=tp, center_1=center_1, center_2=center_2)
    )

    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {outer_1} {center_1} {inner_3}
    .polygon {outer_1} {inner_3} {center_2}
    .polygon {outer_1} {inner_4} {center_1}
    .polygon {outer_1} {center_2} {inner_4}
    .polygon {outer_2} {center_1} {inner_4}
    .polygon {outer_2} {inner_4} {center_2}
    .polygon {outer_2} {inner_5} {center_1}
    .polygon {outer_2} {center_2} {inner_5}
    .polygon {outer_3} {center_1} {inner_5}
    .polygon {outer_3} {inner_5} {center_2}
    .polygon {outer_3} {inner_1} {center_1}
    .polygon {outer_3} {center_2} {inner_1}
    .polygon {outer_4} {center_1} {inner_1}
    .polygon {outer_4} {inner_1} {center_2}
    .polygon {outer_4} {inner_2} {center_1}
    .polygon {outer_4} {center_2} {inner_2}
    .polygon {outer_5} {center_1} {inner_2}
    .polygon {outer_5} {inner_2} {center_2}
    .polygon {outer_5} {inner_3} {center_1}
    .polygon {outer_5} {center_2} {inner_3}
    """.format(
        **color_and_points
    )
    add_vrml_model(bild, "talaia " + message, message)


def pentagon(center, n, p, k, size, color1, rt, name, message, transparency=0):
    """
    Used for residues H and W.
    Describe a pentagon shape placed on top of the side chain's residue.
    Colour depends on residue type.
    """
    n = Point(*n)
    center = Point(*center)
    p = Point(*p)
    k = Point(*k)
    tp = transparency
    color1 = color1[rt]
    vec_cn = n - center
    vec_kp = k - p
    thickness = size / 4.5

    perp_for = normalize(cross(vec_cn, vec_kp)) * thickness
    perp_back = -1 * perp_for
    c2 = (vec_cn) + center
    o1 = (vec_kp) + center
    d1 = c2 - center

    vertex = [
        _rotate(o1, center, c2, alpha, d1) + center
        for alpha in (72, 144, 216, 288, 360)
    ]

    outer = []

    for v in vertex:
        outer.append((-1 * vec_cn) + v)

    ct = (-1 * vec_cn) + center

    color_and_points = dict(
        front_1=perp_for + outer[0],
        front_2=perp_for + outer[1],
        front_3=perp_for + outer[2],
        front_4=perp_for + outer[3],
        front_5=perp_for + outer[4],
        back_1=perp_back + outer[0],
        back_2=perp_back + outer[1],
        back_3=perp_back + outer[2],
        back_4=perp_back + outer[3],
        back_5=perp_back + outer[4],
        center_1=perp_for + ct,
        center_2=perp_back + ct,
        color1=color1,
        tp=tp,
    )

    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {front_1} {front_2} {center_1}
    .polygon {front_3} {center_1} {front_2}
    .polygon {front_1} {center_1} {front_5}
    .polygon {front_3} {front_4} {center_1}
    .polygon {front_5} {center_1} {front_4}
    .polygon {back_1} {center_2} {back_2}
    .polygon {back_3} {back_2} {center_2}
    .polygon {back_1} {back_5} {center_2}
    .polygon {back_3} {center_2} {back_4}
    .polygon {back_5} {back_4} {center_2}
    .polygon {back_1} {back_2} {front_1}
    .polygon {back_2} {front_2} {front_1}
    .polygon {back_2} {back_3} {front_2}
    .polygon {back_4} {front_4} {front_3}
    .polygon {back_3} {front_3} {front_2}
    .polygon {back_3} {back_4} {front_3}
    .polygon {back_4} {back_5} {front_4}
    .polygon {back_5} {front_5} {front_4}
    .polygon {back_5} {back_1} {front_5}
    .polygon {back_1} {front_1} {front_5}
    """.format(
        **color_and_points
    )
    add_vrml_model(bild, "talaia " + message, message)


def hexagon(center, n, p, k, u, size, color1, name, message, transparency=0):
    """
    Used for residues F and Y.
    Describe a hexagon shape placed on the residues ring. Colour depends
    on residue type.
    """
    n = Point(*n)
    center = Point(*center)
    p = Point(*p)
    k = Point(*k)
    u = Point(*u)
    tp = transparency
    if isinstance(color1, tuple):
        color2 = color1[1]
        color1 = color1[0]
    else:
        color1 = color1
        color2 = color1

    vec_kn = k - n
    vec_kp = k - p
    half_ring = center.distance(u) / (-2 * k.distance(p))
    thickness = size / 4.5
    vec_for = normalize(cross(vec_kn, vec_kp)) * thickness
    vec_back = -1 * vec_for
    if size == 3.5:  #TRP
        adjust1 = 1 / length(vec_kn) * 1.7
        adj_vec_kn = adjust1 * vec_kn
        adjust2 = 1 / length(vec_kp) * 1.7
        adj_vec_kp = adjust2 * vec_kp
        c2 = adj_vec_kn + center
        o1 = adj_vec_kp + center
    else:  #PHE and TYR
        c2 = vec_kn + center
        o1 = vec_kp + center
    d1 = center - o1
    vertex = [
        _rotate(o1, center, c2, alpha, d1) + center
        for alpha in (0, 45, 135, 180, 225, 315)
    ]
    outer = []
    for v in vertex:
        outer.append((vec_kp * half_ring) + v)
    color_and_points = dict(
        front_1=vec_for + outer[0],
        front_2=vec_for + outer[1],
        front_3=vec_for + outer[2],
        front_4=vec_for + outer[3],
        front_5=vec_for + outer[4],
        front_6=vec_for + outer[5],
        back_1=vec_back + outer[0],
        back_2=vec_back + outer[1],
        back_3=vec_back + outer[2],
        back_4=vec_back + outer[3],
        back_5=vec_back + outer[4],
        back_6=vec_back + outer[5],
        center_1=vec_for + ((vec_kp * half_ring) + center),
        center_2=vec_back + ((vec_kp * half_ring) + center),
        color1=color1,
        color2=color2,
        tp=tp,
    )

    vec32 = (color_and_points['front_3'] - color_and_points['front_2']) / 3

    mid_1 = vec32 + color_and_points['front_2']
    mid_2 = vec32 + color_and_points['front_6']
    mid_3 = vec32 + color_and_points['back_2']
    mid_4 = vec32 + color_and_points['back_6']

    color_and_points.update(dict(mid_1=mid_1, mid_2=mid_2, mid_3=mid_3, mid_4=mid_4))

    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {mid_1} {mid_2} {front_3}
    .polygon {mid_2} {front_3} {front_5}
    .polygon {front_3} {front_4} {front_5}
    .polygon {mid_3} {mid_4} {back_3}
    .polygon {mid_4} {back_3} {back_5}
    .polygon {back_3} {back_4} {back_5}
    .polygon {mid_1} {mid_3} {back_3}
    .polygon {mid_1} {back_3} {front_3}
    .polygon {mid_3} {mid_4} {back_5}
    .polygon {mid_2} {back_5} {front_5}
    .polygon {mid_2} {mid_4} {back_5}
    .polygon {front_3} {front_4} {back_3}
    .polygon {back_3} {back_4} {front_4}
    .polygon {front_4} {front_5} {back_5}
    .polygon {back_4} {back_5} {front_4}
    .color {color2}
    .polygon {front_1} {front_2} {front_6}
    .polygon {back_1} {back_2} {back_6}
    .polygon {front_6} {front_2} {mid_1}
    .polygon {front_6} {mid_1} {mid_2}
    .polygon {back_6} {back_2} {mid_3}
    .polygon {back_6} {mid_3} {mid_4}
    .polygon {front_1} {front_6} {back_1}
    .polygon {front_6} {back_6} {back_1}
    .polygon {front_1} {front_2} {back_1}
    .polygon {front_2} {back_2} {back_1}
    .polygon {front_6} {mid_2} {mid_4}
    .polygon {front_6} {back_6} {mid_4}
    .polygon {front_2} {mid_1} {mid_3}
    .polygon {front_2} {back_2} {mid_3}
    """.format(
        **color_and_points
    )
    add_vrml_model(bild, "talaia " + message, message)


def add_vrml_model(vrml_string, name, message, prefix_msg="#? "):
    """
    create vrml model from bild string
    """
    f = StringIO(dedent(vrml_string))
    try:
        vrml = openBildFileObject(f, "<string>", "vrml_" + name)
        vrml[0].molecule = BalloonPopup(prefix_msg + message)
    except chimera.NotABug:
        print vrml_string
    else:
        chimera.openModels.add(vrml)
        return vrml

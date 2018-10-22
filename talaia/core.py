#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
TALAIA
======

Simplistic 3D dictionary of amino acids using geometric shapes
for UCSF Chimera.
"""


from textwrap import dedent
import math
import argparse
from collections import defaultdict

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

import numpy as np

import chimera
from Bld2VRML import openFileObject as openBildFileObject
from chimera import runCommand as run, cross, Point, Vector, preferences
from chimera import selection, specifier
from Shape import shapecmd
import Matrix as M


AMINOACIDS = {
    "GLY": dict(
        name="GLY",
        figure="sphere",
        size=0.9,
        color1=(0, 1, 0),
        all_atoms=["N", "O", "CA", "C"],
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
        color1="orange",
        all_atoms=["N", "CA", "O", "C", "CB", "OG"],
        center="OG",
        n="CA",
        p="O",
    ),
    "THR": dict(
        name="THR",
        figure="cone",
        size=3.4,
        color1="orange",
        all_atoms=["N", "CA", "O", "C", "CB", "OG1", "CG2"],
        center="OG1",
        n="CA",
        p="O",
    ),
    "CYS": dict(
        name="CYS",
        figure="cone",
        size=3.4,
        color1="orange",
        all_atoms=["N", "CA", "C", "O", "CB", "SG"],
        center="SG",
        n="CA",
        p="O",
    ),
    "ASN": dict(
        name="ASN",
        figure="triangle",
        size=3.4,
        color1="orange",
        color2="blue",
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
        k="CB",
        n="CG",
        p="OD1",
        u="ND2",
    ),
    "GLN": dict(
        name="GLN",
        figure="triangle",
        size=3.4,
        color1="orange",
        color2="blue",
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
        k="CA",
        n="CD",
        p="OE1",
        u="NE2",
    ),
    "TYR": dict(
        name="TYR",
        figure="hexagon",
        size=3.4,
        color1="orange",
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
        color2="dark red",
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
        color2="dark red",
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
        color1="yellow",
        all_atoms=["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
        center="CG",
        n="CB",
        k="CD2",
        p="ND1",
        rot_factor=-1,
    ),
}

ALTERNATIVE_NAME = {"CYX": "CYS", "HIE": "HIS", "HE1": "HIS", "CS1": "CYS"}

IONS_AND_METALS = ["MG", "NI", "PT1", "K", "CA", "NA", "F"]


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
            or chimera.specifier.evalSpec("all").residues()
        )

    elif isinstance(selection, basestring):
        nearRes = chimera.specifier.evalSpec(selection).residues()
    else:  # we assume it is a list of residues
        nearRes = selection
    for res in nearRes:
        res_info = AMINOACIDS.get(res.type)
        if res_info is None:
            if res.type in ALTERNATIVE_NAME.keys():
                alt_name = ALTERNATIVE_NAME.get(res.type)
                res_info = AMINOACIDS.get(alt_name)
            elif res.isMetal:
                for at in res.atoms:
                    ion_center = at.coord()
                    ion_star(ion_center, 3 * at.radius, at.color)
                    continue
            elif res.isHet:
                continue
            else:
                if "CA" not in res.atomsMap:
                    center = res.atoms[0].coord()
                    continue

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
            rt = res_info.get("rot_factor")

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
                shapecmd.sphere_shape(radius=r1, center=q, color=(0, 1, 0, 1 - tp))
            elif shape == "ellipsoid":
                ellipsoid(k, n, tp)
            elif shape == "triangle":
                triangle(n, p, k, u, r1, color1, color2, tp, name)
            elif shape == "cone":
                cone(center, n, p, r1, color1, name, tp)
            elif shape == "cube":
                cube(center, n, p, r1, color1, name, tp)
            elif shape == "rectangle":
                rectangle(center, n, p, k, r1, color1, color2, name, tp)
            elif shape == "pentagon":
                pentagon(center, n, p, k, r1, color1, rt, name, tp)
            elif shape == "hexagon":
                hexagon(center, n, p, k, u, r1, color1, name, tp)

    subscribe_events()


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
        ):
            chimera.openModels.close([item])
        if hasattr(item, "_depicted"):
            del item._depicted
            for r in item.residues:
                if hasattr(r, "_depicted"):
                    del r._depicted

    unsubscribe_events()


def subscribe_events():
    chimera._depicter_handler = chimera.triggers.addHandler("Molecule", callback, None)


def unsubscribe_events():
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
            deleted_mols,
        ]
    ):
        disable()
        if options:
            enable(**options)


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


def ellipsoid(k, n, transparency=0):
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
    shapecmd.sphere_shape(
        radius=(2, 1, 1), center=q, rotation=(x, y, z, theta), color=(0, 1, 0, tp)
    )


def cube(center, n, p, size, color1, name, transparency=0):
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
    add_vrml_model(bild, "cube_" + name)


def triangle(n, p, k, u, size, color1, color2, transparency, name):
    """
    Used for residues E, D, Q, N.
    Define oriented orange triangle with one side of different
    color depending on residue type.
    """
    n = Point(*n)
    p = Point(*p)
    k = Point(*k)
    u = Point(*u)
    tp = transparency
    # create two vectors and adapt their length to calculate center point.
    vec_nk = n - k
    vec_pu = p - u
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
    # o4 = perp2 + x2

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

    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {s1} {s2} {s3}
    .polygon {s2} {s3} {s4}
    .polygon {s1} {s2} {s5}
    .polygon {s2} {s6} {s5}
    .color {color2}
    .polygon {s6} {s5} {d1}
    .polygon {s6} {d1} {d2}
    .color {color2}
    .polygon {s3} {d2} {d1}
    .polygon {s3} {s4} {d2}
    .color {color1}
    .polygon {s1} {s3} {s5}
    .polygon {s2} {s4} {s6}
    """.format(
        color1=color1, color2=color2, tp=tp, **points
    )
    add_vrml_model(bild, "triangle_" + name)


def cone(center, n, p, size, color1, name, transparency=0):
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
    color_and_points.update(dict(x1=x1, x4=x4, color1=color1, tp=tp))

    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {x4} {outer_1} {outer_2}
    .polygon {x4} {outer_2} {outer_3}
    .polygon {x4} {outer_3} {outer_1}
    .polygon {outer_1} {outer_2} {outer_3}
    """.format(
        **color_and_points
    )
    add_vrml_model(bild, "cone_" + name)


def rectangle(center, n, p, k, size, color1, color2, name, transparency=0):
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
    add_vrml_model(bild, "rectangle_" + name)


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


def star(center, n, p, size, color1, transparency=0):
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
    add_vrml_model(bild, "star")


def pentagon(center, n, p, k, size, color1, rt, name, transparency=0):
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

    vec_cn = n - center
    vec_kp = k - p
    # half_length = size / 2.0
    thickness = size / 4.5

    perp_for = normalize(cross(vec_cn, vec_kp)) * thickness
    perp_back = -1 * perp_for
    c2 = (vec_cn * abs(1 / rt)) + center
    o1 = (vec_kp * abs(1 / rt)) + center
    d1 = c2 - center

    vertex = [
        _rotate(o1, center, c2, alpha, d1) + center
        for alpha in (72, 144, 216, 288, 360)
    ]

    outer = []

    for v in vertex:
        outer.append(((rt / abs(rt)) * vec_cn) + v)

    ct = (rt / abs(rt)) * vec_cn + center

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
    add_vrml_model(bild, "pentagon_" + name)


def hexagon(center, n, p, k, u, size, color1, name, transparency=0):
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
    vec_kn = k - n
    vec_kp = k - p
    half_ring = center.distance(u) / (-2 * k.distance(p))
    thickness = size / 4.5
    vec_for = normalize(cross(vec_kn, vec_kp)) * thickness
    vec_back = -1 * vec_for
    if size == 3.5:
        adjust1 = 1 / length(vec_kn) * 1.7
        adj_vec_kn = adjust1 * vec_kn
        adjust2 = 1 / length(vec_kp) * 1.7
        adj_vec_kp = adjust2 * vec_kp
        c2 = adj_vec_kn + center
        o1 = adj_vec_kp + center
    else:
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
        tp=tp,
    )
    bild = """
    .color {color1}
    .transparency {tp}
    .polygon {front_1} {front_2} {center_1}
    .polygon {front_1} {center_1} {front_6}
    .polygon {front_3} {center_1} {front_2}
    .polygon {front_3} {front_4} {center_1}
    .polygon {front_5} {center_1} {front_4}
    .polygon {front_5} {front_6} {center_1}
    .polygon {back_1} {center_2} {back_2}
    .polygon {back_1} {back_6} {center_2}
    .polygon {back_3} {back_2} {center_2}
    .polygon {back_3} {center_2} {back_4}
    .polygon {back_5} {back_4} {center_2}
    .polygon {back_5} {center_2} {back_6}
    .polygon {back_1} {back_2} {front_1}
    .polygon {back_2} {front_2} {front_1}
    .polygon {back_2} {back_3} {front_2}
    .polygon {back_3} {front_3} {front_2}
    .polygon {back_3} {back_4} {front_3}
    .polygon {back_4} {front_4} {front_3}
    .polygon {back_4} {back_5} {front_4}
    .polygon {back_5} {front_5} {front_4}
    .polygon {back_5} {back_6} {front_5}
    .polygon {back_6} {front_6} {front_5}
    .polygon {back_6} {back_1} {front_6}
    .polygon {back_1} {front_1} {front_6}
    """.format(
        **color_and_points
    )
    add_vrml_model(bild, "hexagon_" + name)


def add_vrml_model(vrml_string, name):
    f = StringIO(dedent(vrml_string))
    try:
        vrml = openBildFileObject(f, "<string>", "vrml_" + name)
    except chimera.NotABug:
        print(vrml_string)
    else:
        chimera.openModels.add(vrml)
        return vrml

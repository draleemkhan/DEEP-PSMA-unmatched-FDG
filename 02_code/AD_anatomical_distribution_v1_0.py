"""
AD_anatomical_distribution_v1_0.py - PAPER 3, SAP section 8.3: anatomical distribution of the primary unmatched FDG
burden, and the anatomical consistency flag of matched pairs (separate descriptive audit).

Implements: locked consolidated Step 2 SAP (MD5 0c64711db1207e31a70cf3f558ccbd7d) section 8.3; dated implementation
clarification 00_protocol\\StepSAP_anatomical_distribution_implementation_clarification_20261112.md
(MD5 18c59df4c97b7cdba9fabb45fc458fb0; author decisions A1-A8 and two safeguards). Locked results are read only.

Distribution units: the 217 primary V_U components (P3_conn26_20261112_001915). Native FDG CT TotalSegmentator labels
(totseg_24); native voxel centres (P1 label image, 26-connectivity) -> nearest FDG CT voxel by floor(x + 0.5).
Assignment (strictly > 50 % of ALL native voxel centres), in order:
  inside CT <= 0.5 -> site-indeterminate: anatomy unavailable
  bone > 0.5 -> bone: spinal > 0.5 -> spinal region; else one other skeletal group > 0.5; else mixed skeletal groups
  one visceral organ > 0.5 -> that organ
  non-organ soft tissue > 0.5 -> soft tissue (location, never "nodal")
  otherwise -> site-indeterminate: anatomically ambiguous
Spinal regions: LOCKED PSMA-space plane table applied to the T_A-mapped mean z of the native voxel centres
  (TA_centroid_z of the primary P2 run = the R6 convention): cervical z >= C7/T1; thoracic T12/L1 <= z < C7/T1;
  lumbar L5/S1 <= z < T12/L1; sacral z < L5/S1; on a plane -> superior region; UNDEFINED plane -> region unassigned.
Soft-tissue location: T_A-mapped z > C2/C3 -> unassignable (head); T12/L1 undefined -> location unassigned;
  z >= T12/L1 -> supradiaphragmatic location (operational proxy, not the anatomical diaphragm); below: native mean z
  (P1 centroid) vs the inferior end of the native FDG CT aorta label (bifurcation): >= -> retroperitoneal level;
  < -> iliac zone (region unassigned), unless below the inferior end of the pelvic bones -> unassignable (below
  pelvis); aorta absent -> infradiaphragmatic level unassigned; pelvic bones absent -> no lower limit (recorded).
Consistency flag (separate audit): primary V_M components x listed matched PSMA components; coarse compartment (bone /
  visceral organ / soft tissue / not evaluable) of each side by the same rule (PSMA side: PSMA TTB components,
  26-connectivity, raster order of first voxel = MATLAB bwconncomp; native PSMA CT labels). Flag = both evaluable and
  different. ID check: every matched FDG component with final distance 0 shares >= 1 voxel with a listed PSMA component.
Outputs: 06_output\\AD_[PREFLIGHT_]<run>\\ AD_components, AD_categories, AD_flag_pairs, AD_flag_summary (CSV),
  AD_summary (txt), AD_manifest (json).
Usage (folder 02_code, env totalseg):
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u AD_anatomical_distribution_v1_0.py --selftest
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u AD_anatomical_distribution_v1_0.py --preflight
  & "C:\\Users\\Administrator\\miniconda3\\envs\\totalseg\\python.exe" -u AD_anatomical_distribution_v1_0.py
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

DATA_ROOT = Path(r"D:\DEEP PSMA\DEEP PSMA EXTRACTED")
OUT_DIR = Path(r"D:\DEEP PSMA\PAPER3_Discordance\06_output")
PROTO_DIR = OUT_DIR.parent / "00_protocol"
P1_RUN, P2_RUN, P2B_RUN, P3_RUN = ("P1_20261111_233936", "P2_conn26_20261111_235603", "P2b_conn26_20261112_000903",
                                   "P3_conn26_20261112_001915")
CONN = 26
INPUTS = {
    (OUT_DIR, f"{P1_RUN}/P1_components_20261111_233936.csv"): "ffafbc7917e1b595eadacafab247fe88",
    (OUT_DIR, f"{P1_RUN}/P1_patients_20261111_233936.csv"): "f9b5b5908e523d1122547ef07d5c0afd",
    (OUT_DIR, f"{P2_RUN}/P2_components_20261111_235603.csv"): "9a649e60c49dafde2e4b3d1b1e35fa3c",
    (OUT_DIR, f"{P2B_RUN}/P2b_psma_pet_grid_20261112_000903.csv"): "ad5504fe4fb4118b727c167b53019983",
    (OUT_DIR, f"{P3_RUN}/P3_components_20261112_001915.csv"): "ac0048cfe059c714b3c293fbe7e88279",
    (OUT_DIR, "SAP_R1_plane_table_20261111_214952.csv"): "db6fe2142e5c6fe7fe8ee6d8f9bc0f54",
    (OUT_DIR, "Audit1B_label_level_20261029_134401.csv"): "896a1681e90da4fe793a3a449d737fde",
    (PROTO_DIR, "StepSAP_consolidated_LOCKED_20261111.md"): "0c64711db1207e31a70cf3f558ccbd7d",
    (PROTO_DIR, "StepSAP_anatomical_distribution_implementation_clarification_20261112.md"): "18c59df4c97b7cdba9fabb45fc458fb0",
}
SPINAL = ["sacrum"] + [f"vertebrae_{s}" for s in ["S1", "L5", "L4", "L3", "L2", "L1"] + [f"T{i}" for i in range(12, 0, -1)]
                                                    + [f"C{i}" for i in range(7, 0, -1)]]
SKELETAL = {"skull": ["skull"],
            "ribs and costal cartilages": [f"rib_{s}_{i}" for s in ("left", "right") for i in range(1, 13)] + ["costal_cartilages"],
            "sternum": ["sternum"], "clavicles": ["clavicula_left", "clavicula_right"],
            "scapulae": ["scapula_left", "scapula_right"], "humeri": ["humerus_left", "humerus_right"],
            "pelvic bones": ["hip_left", "hip_right"], "femurs": ["femur_left", "femur_right"]}
ORGANS = {"spleen": ["spleen"], "kidneys": ["kidney_right", "kidney_left", "kidney_cyst_left", "kidney_cyst_right"],
          "gallbladder": ["gallbladder"], "liver": ["liver"], "stomach": ["stomach"], "pancreas": ["pancreas"],
          "adrenal glands": ["adrenal_gland_right", "adrenal_gland_left"],
          "lungs": ["lung_upper_lobe_left", "lung_lower_lobe_left", "lung_upper_lobe_right", "lung_middle_lobe_right",
                    "lung_lower_lobe_right"],
          "oesophagus": ["esophagus"], "trachea": ["trachea"], "thyroid": ["thyroid_gland"], "small bowel": ["small_bowel"],
          "duodenum": ["duodenum"], "colon": ["colon"], "urinary bladder": ["urinary_bladder"], "prostate": ["prostate"],
          "heart": ["heart"], "brain": ["brain"]}
COARSE = {"bone": "bone", "organ": "visceral organ", "soft": "soft tissue"}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def build_ids(name_to_id):
    """label-id sets from structure names (KeyError if a name is missing from the dictionary)."""
    g = dict(spinal={name_to_id[n] for n in SPINAL},
             skeletal={k: {name_to_id[n] for n in v} for k, v in SKELETAL.items()},
             organs={k: {name_to_id[n] for n in v} for k, v in ORGANS.items()})
    g["bone"] = g["spinal"] | set().union(*g["skeletal"].values())
    g["organ_all"] = set().union(*g["organs"].values())
    return g


# ------------------------------ pure helpers (tested) ------------------------------
def ct_index(pts, origin, spacing, direction, size):
    """nearest CT voxel (k, j, i) by floor(x + 0.5); rows outside the grid -> inside False (as ST2a)."""
    D = np.asarray(direction, float).reshape(3, 3)
    ci = np.linalg.solve(D * np.asarray(spacing, float), (np.atleast_2d(pts) - np.asarray(origin, float)).T).T
    ijk = np.floor(ci + 0.5).astype(np.int64)
    inside = np.all((ijk >= 0) & (ijk < np.asarray(size, np.int64)), axis=1)
    return ijk[:, ::-1], inside


def labels_at(arr, kji, inside):
    out = np.zeros(len(kji), np.int64)
    if inside.any():
        k = kji[inside]
        out[inside] = arr[k[:, 0], k[:, 1], k[:, 2]]
    return out


def tissue_class(lab, inside, G):
    """(class, sub): class in anatomy unavailable / bone / organ / soft / ambiguous. Strictly > 50 % of ALL points."""
    n = len(lab)
    if inside.sum() / n <= 0.5:
        return "unavailable", ""
    li = lab[inside]
    nb = np.isin(li, list(G["bone"])).sum()
    if nb / n > 0.5:
        if np.isin(li, list(G["spinal"])).sum() / n > 0.5:
            return "bone", "spinal"
        for k, ids in G["skeletal"].items():
            if np.isin(li, list(ids)).sum() / n > 0.5:
                return "bone", k
        return "bone", "mixed skeletal groups"
    for k, ids in G["organs"].items():
        if np.isin(li, list(ids)).sum() / n > 0.5:
            return "organ", k
    ns = len(li) - nb - np.isin(li, list(G["organ_all"])).sum()
    if ns / n > 0.5:
        return "soft", ""
    return "ambiguous", ""


def spinal_region(z, P):
    if z >= P["z_C7T1"]:
        return "cervical"
    if not np.isfinite(P["z_T12L1"]):
        return "region unassigned (plane undefined)"
    if z >= P["z_T12L1"]:
        return "thoracic"
    if not np.isfinite(P["z_pelvis_L5S1"]):
        return "region unassigned (plane undefined)"
    if z >= P["z_pelvis_L5S1"]:
        return "lumbar"
    return "sacral"


def soft_location(z_ta, z_nat, P, bif_z, pelvis_min_z):
    if z_ta > P["z_head_C2C3"]:
        return "unassignable soft tissue", "head"
    if not np.isfinite(P["z_T12L1"]):
        return "soft tissue", "location unassigned (plane undefined)"
    if z_ta >= P["z_T12L1"]:
        return "soft tissue", "supradiaphragmatic location"
    if not np.isfinite(bif_z):
        return "soft tissue", "infradiaphragmatic level unassigned"
    if z_nat >= bif_z:
        return "soft tissue", "retroperitoneal level"
    if np.isfinite(pelvis_min_z) and z_nat < pelvis_min_z:
        return "unassignable soft tissue", "below pelvis"
    return "soft tissue", "iliac zone (region unassigned)"


def category(cls, sub, z_ta, z_nat, P, bif_z, pelvis_min_z):
    if cls == "unavailable":
        return "site-indeterminate", "anatomy unavailable"
    if cls == "ambiguous":
        return "site-indeterminate", "anatomically ambiguous"
    if cls == "bone":
        return ("bone", "spinal: " + spinal_region(z_ta, P)) if sub == "spinal" else ("bone", sub)
    if cls == "organ":
        return "visceral organ", sub
    return soft_location(z_ta, z_nat, P, bif_z, pelvis_min_z)


def compartment(cls):
    return COARSE.get(cls, "not evaluable")


def label_first_voxel_order(bw):
    """26-connected components numbered by the raster (C) order of their first voxel."""
    lab, n = ndimage.label(bw, structure=np.ones((3, 3, 3), bool))
    if n == 0:
        return lab.astype(np.int32), 0
    flat = lab.ravel()
    ids, first = np.unique(flat, return_index=True)
    keep = ids > 0
    ids, first = ids[keep], first[keep]
    remap = np.zeros(n + 1, np.int32)
    remap[ids[np.argsort(first)]] = np.arange(1, n + 1, dtype=np.int32)
    return remap[lab], int(n)


def parse_ids(s):
    if s is None or (isinstance(s, float) and np.isnan(s)) or str(s).strip() == "":
        return []
    return [int(float(x)) for x in str(s).split(";") if x.strip()]


# ------------------------------ self-test ------------------------------
def selftest():
    res = []

    def ok(name, cond):
        res.append((name, bool(cond))); print(("PASS " if cond else "FAIL ") + name)

    names = (["spleen", "kidney_right", "kidney_left", "gallbladder", "liver", "stomach", "pancreas", "adrenal_gland_right",
              "adrenal_gland_left", "lung_upper_lobe_left", "lung_lower_lobe_left", "lung_upper_lobe_right",
              "lung_middle_lobe_right", "lung_lower_lobe_right", "esophagus", "trachea", "thyroid_gland", "small_bowel",
              "duodenum", "colon", "urinary_bladder", "prostate", "kidney_cyst_left", "kidney_cyst_right"] + SPINAL
             + ["heart", "aorta", "skull", "brain", "sternum", "costal_cartilages", "clavicula_left", "clavicula_right",
                "scapula_left", "scapula_right", "humerus_left", "humerus_right", "hip_left", "hip_right", "femur_left",
                "femur_right", "iliopsoas_left"] + [f"rib_{s}_{i}" for s in ("left", "right") for i in range(1, 13)])
    nid = {n: i + 1 for i, n in enumerate(dict.fromkeys(names))}
    G = build_ids(nid)
    L = lambda *ns: np.array([nid[x] if isinstance(x, str) else x for x in ns])
    T = lambda n: np.ones(n, bool)
    v = nid["vertebrae_L3"]; r = nid["rib_left_4"]; li = nid["liver"]; ps = nid["iliopsoas_left"]
    ok("outside CT 50 % exactly -> anatomy unavailable", tissue_class(L(v, v, v, v), np.array([1, 1, 0, 0], bool), G)[0] == "unavailable")
    ok("outside CT 25 % -> evaluated", tissue_class(L(v, v, v, v), np.array([1, 1, 1, 0], bool), G) == ("bone", "spinal"))
    ok("bone exactly 50 % -> not bone (strict >)", tissue_class(L(v, v, li, 0), T(4), G)[0] != "bone")
    ok("bone 3/5 spinal 3/5 -> spinal", tissue_class(L(v, v, v, 0, 0), T(5), G) == ("bone", "spinal"))
    ok("bone 4/5, spinal 2/5, ribs 2/5 -> mixed skeletal groups", tissue_class(L(v, v, r, r, 0), T(5), G) == ("bone", "mixed skeletal groups"))
    ok("ribs 3/5 -> ribs and costal cartilages", tissue_class(L(r, r, nid["costal_cartilages"], 0, 0), T(5), G) == ("bone", "ribs and costal cartilages"))
    ok("liver 3/5 -> organ liver", tissue_class(L(li, li, li, 0, v), T(5), G) == ("organ", "liver"))
    ok("lung lobes combined -> lungs", tissue_class(L(nid["lung_upper_lobe_left"], nid["lung_lower_lobe_right"], nid["lung_middle_lobe_right"], 0), T(4), G) == ("organ", "lungs"))
    ok("muscle and unlabelled 3/5 -> soft tissue", tissue_class(L(ps, 0, 0, li, v), T(5), G) == ("soft", ""))
    ok("no class > 50 % -> anatomically ambiguous", tissue_class(L(v, v, li, li, 0, 0), T(6), G) == ("ambiguous", ""))
    ok("outside points count in the denominator (2 bone of 5, 2 outside)", tissue_class(L(v, v, 0, 0, 0), np.array([1, 1, 1, 0, 0], bool), G) == ("ambiguous", ""))
    P = dict(z_head_C2C3=100.0, z_C7T1=40.0, z_T12L1=-200.0, z_pelvis_L5S1=-400.0)
    ok("spine: on C7/T1 -> cervical; just below -> thoracic", spinal_region(40.0, P) == "cervical" and spinal_region(39.9, P) == "thoracic")
    ok("spine: on T12/L1 -> thoracic; on L5/S1 -> lumbar; below -> sacral",
       spinal_region(-200.0, P) == "thoracic" and spinal_region(-400.0, P) == "lumbar" and spinal_region(-400.1, P) == "sacral")
    ok("spine: above C2/C3 still cervical (C1-C2)", spinal_region(150.0, P) == "cervical")
    PU = dict(P, z_T12L1=np.nan); PL = dict(P, z_pelvis_L5S1=np.nan)
    ok("spine: T12/L1 undefined -> unassigned below C7/T1, cervical still assigned",
       spinal_region(0.0, PU).startswith("region unassigned") and spinal_region(50.0, PU) == "cervical")
    ok("spine: L5/S1 undefined -> lumbar/sacral unassigned, thoracic assigned",
       spinal_region(-300.0, PL).startswith("region unassigned") and spinal_region(-100.0, PL) == "thoracic")
    ok("soft: above C2/C3 -> unassignable (head); on C2/C3 -> supradiaphragmatic",
       soft_location(100.1, 0, P, -300, -500) == ("unassignable soft tissue", "head")
       and soft_location(100.0, 0, P, -300, -500) == ("soft tissue", "supradiaphragmatic location"))
    ok("soft: on T12/L1 -> supradiaphragmatic; below -> by native z vs bifurcation",
       soft_location(-200.0, -999, P, -300, -500)[1] == "supradiaphragmatic location"
       and soft_location(-201, -300.0, P, -300, -500)[1] == "retroperitoneal level"
       and soft_location(-201, -300.1, P, -300, -500)[1] == "iliac zone (region unassigned)")
    ok("soft: below pelvic bones -> unassignable (below pelvis); pelvic bones absent -> iliac zone",
       soft_location(-201, -600, P, -300, -500) == ("unassignable soft tissue", "below pelvis")
       and soft_location(-201, -600, P, -300, np.nan)[1] == "iliac zone (region unassigned)")
    ok("soft: aorta absent -> infradiaphragmatic level unassigned; T12/L1 undefined -> location unassigned",
       soft_location(-201, -300, P, np.nan, -500)[1] == "infradiaphragmatic level unassigned"
       and soft_location(0, 0, PU, -300, -500)[1] == "location unassigned (plane undefined)")
    ok("soft tissue is never labelled nodal", all("nodal" not in "".join(soft_location(z, zn, P, -300, -500))
                                                 for z in (150, 0, -250) for zn in (-100, -350, -600)))
    ok("category: site-indeterminate kinds kept separate",
       category("unavailable", "", 0, 0, P, -300, -500) == ("site-indeterminate", "anatomy unavailable")
       and category("ambiguous", "", 0, 0, P, -300, -500) == ("site-indeterminate", "anatomically ambiguous"))
    ok("category: spinal region from T_A z; organ; skeletal group",
       category("bone", "spinal", -300, 999, P, -300, -500) == ("bone", "spinal: lumbar")
       and category("organ", "liver", 0, 0, P, 0, 0) == ("visceral organ", "liver")
       and category("bone", "pelvic bones", 0, 0, P, 0, 0) == ("bone", "pelvic bones"))
    ok("compartments: bone / visceral organ / soft tissue / not evaluable",
       [compartment(c) for c in ("bone", "organ", "soft", "unavailable", "ambiguous")]
       == ["bone", "visceral organ", "soft tissue", "not evaluable", "not evaluable"])
    bw = np.zeros((4, 5, 6), bool); bw[3, 4, 5] = True; bw[0, 0, 2] = True; bw[0, 0, 3] = True; bw[2, 2, 0] = True
    lab, n = label_first_voxel_order(bw)
    ok("PSMA components numbered by raster order of first voxel (26-conn)", n == 3 and lab[0, 0, 2] == 1 and lab[0, 0, 3] == 1
       and lab[2, 2, 0] == 2 and lab[3, 4, 5] == 3)
    bw2 = np.zeros((3, 3, 3), bool); bw2[0, 0, 0] = True; bw2[1, 1, 1] = True
    ok("26-connectivity joins diagonal voxels", label_first_voxel_order(bw2)[1] == 1)
    ok("matched-list parsing", parse_ids("47;51") == [47, 51] and parse_ids(4) == [4] and parse_ids(np.nan) == [] and parse_ids("") == [])
    kji, ins = ct_index(np.array([[0.4, 0.0, 0.0], [0.6, 0.0, 0.0], [-0.6, 0, 0]]), (0, 0, 0), (1, 1, 1), (1, 0, 0, 0, 1, 0, 0, 0, 1), (3, 3, 3))
    ok("CT nearest voxel by floor(x + 0.5); outside detected", kji[0, 2] == 0 and kji[1, 2] == 1 and not ins[2] and ins[0])
    n_ok = sum(c for _, c in res)
    print(f"\nSELF-TEST {n_ok}/{len(res)} passed")
    return n_ok == len(res)


# ------------------------------ run ------------------------------
def phys(img, kji):
    D = np.array(img.GetDirection()).reshape(3, 3); S = np.array(img.GetSpacing())
    return np.array(img.GetOrigin()) + (np.asarray(kji)[:, ::-1].astype(float) * S) @ D.T


def label_extent_min_z(img, arr, ids):
    kji = np.argwhere(np.isin(arr, list(ids)))
    return float(phys(img, kji)[:, 2].min()) if len(kji) else np.nan


def process(cid, VU, VM, G, nid, P1C, P1P, P2C, PL, GR, imgs):
    import SimpleITK as sitk
    dF, dP = DATA_ROOT / cid / "FDG", DATA_ROOT / cid / "PSMA"
    lab_path = OUT_DIR / P1_RUN / f"P1_{cid}_FDG_conn{CONN}.nii.gz"
    if md5(lab_path) != P1P.loc[cid, "label_image_md5"]:
        raise RuntimeError(f"{cid}: P1 label image MD5 differs from the P1 patient table")
    lab_img = sitk.ReadImage(str(lab_path)); lab = sitk.GetArrayFromImage(lab_img)
    idx = ndimage.value_indices(lab, ignore_value=0)
    ctF = sitk.ReadImage(str(dF / "totseg_24.nii.gz")); aF = sitk.GetArrayFromImage(ctF).astype(np.int64)
    imgs[f"{cid}/FDG/totseg_24.nii.gz"] = md5(dF / "totseg_24.nii.gz")
    gF = (ctF.GetOrigin(), ctF.GetSpacing(), ctF.GetDirection(), ctF.GetSize())
    bif = label_extent_min_z(ctF, aF, {nid["aorta"]})
    pelv = label_extent_min_z(ctF, aF, {nid["hip_left"], nid["hip_right"]})
    P = PL.loc[cid]; Pd = {k: (float(P[k]) if pd.notna(P[k]) else np.nan) for k in ("z_head_C2C3", "z_C7T1", "z_T12L1", "z_pelvis_L5S1")}

    def fdg_class(q):
        kji = np.column_stack(idx[q])
        if len(kji) != int(P1C.loc[(cid, q), "voxels"]):
            raise RuntimeError(f"{cid} c{q}: native voxel count differs from P1")
        k2, ins = ct_index(phys(lab_img, kji), *gF)
        return tissue_class(labels_at(aF, k2, ins), ins, G), float(ins.mean())

    rows, pairs = [], []
    for r in VU.itertuples():
        q = int(r.component_id)
        (cls, sub), fin = fdg_class(q)
        zta = float(P2C.loc[(cid, q), "TA_centroid_z"]); znat = float(P1C.loc[(cid, q), "centroid_z_LPS_mm"])
        cat, sc = category(cls, sub, zta, znat, Pd, bif, pelv)
        rows.append(dict(case_id=cid, component_id=q, volume_native_mL=float(r.volume_native_mL), frac_inside_CT=fin,
                         tissue_class=cls, tissue_sub=sub, category=cat, subcategory=sc, TA_mean_z=zta, native_mean_z=znat,
                         z_C2C3=Pd["z_head_C2C3"], z_C7T1=Pd["z_C7T1"], z_T12L1=Pd["z_T12L1"], z_L5S1=Pd["z_pelvis_L5S1"],
                         aorta_inferior_z=bif, pelvic_bone_inferior_z=pelv))
    if len(VM):
        ttb = sitk.ReadImage(str(dP / "TTB.nii.gz")); ctP = sitk.ReadImage(str(dP / "totseg_24.nii.gz"))
        imgs[f"{cid}/PSMA/TTB.nii.gz"] = md5(dP / "TTB.nii.gz"); imgs[f"{cid}/PSMA/totseg_24.nii.gz"] = md5(dP / "totseg_24.nii.gz")
        aP = sitk.GetArrayFromImage(ctP).astype(np.int64); gP = (ctP.GetOrigin(), ctP.GetSpacing(), ctP.GetDirection(), ctP.GetSize())
        plab, npc = label_first_voxel_order(sitk.GetArrayFromImage(ttb) > 0)
        pidx = ndimage.value_indices(plab, ignore_value=0)
        vs = None
        if (VM.d_final_mm == 0).any():
            vp = OUT_DIR / P2B_RUN / f"P2b_voxelsets_{cid}.csv"
            if md5(vp) != GR.loc[cid, "csv_md5"]:
                raise RuntimeError(f"{cid}: voxel-set CSV MD5 differs from the P2b grid table")
            vs = pd.read_csv(vp); vs = vs[vs.role == 0]
        pcache = {}
        for r in VM.itertuples():
            q = int(r.component_id); ids = parse_ids(r.psma_components_matched)
            if not ids or max(ids) > npc:
                raise RuntimeError(f"{cid} c{q}: matched PSMA component list invalid ({r.psma_components_matched}; {npc} components)")
            if r.d_final_mm == 0:
                lin = vs.loc[vs.component_id == q, "linear_index_1based"].to_numpy() - 1
                if not np.isin(plab.ravel()[lin], ids).any():
                    raise RuntimeError(f"{cid} c{q}: d = 0 but no listed PSMA component shares a voxel (ID convention)")
            (cf, sf), _ = fdg_class(q)
            for p in ids:
                if p not in pcache:
                    kp = np.column_stack(pidx[p]); k2, ins = ct_index(phys(ttb, kp), *gP)
                    pcache[p] = tissue_class(labels_at(aP, k2, ins), ins, G)
                cp, sp = pcache[p]
                a, b = compartment(cf), compartment(cp)
                pairs.append(dict(case_id=cid, fdg_component=q, psma_component=p, volume_native_mL=float(r.volume_native_mL),
                                  fdg_compartment=a, fdg_detail=sf, psma_compartment=b, psma_detail=sp,
                                  evaluable=(a != "not evaluable" and b != "not evaluable"),
                                  flagged=(a != "not evaluable" and b != "not evaluable" and a != b)))
    return rows, pairs


def cell(d):
    return dict(components=int(len(d)), native_mL=float(d.volume_native_mL.sum()), patients=int(d.case_id.nunique()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); ap.add_argument("--preflight", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not selftest():
        sys.exit("self-test failed: no run")
    t0 = time.time()
    bad = [f"{d / n}: {md5(d / n)} != {m}" for (d, n), m in INPUTS.items() if md5(d / n) != m]
    if bad:
        sys.exit("INPUT MD5 MISMATCH:\n" + "\n".join(bad))
    A1B = pd.read_csv(OUT_DIR / "Audit1B_label_level_20261029_134401.csv")
    dic = A1B[A1B.tracer == "FDG"][["label_id", "structure"]].drop_duplicates()
    if len(dic) != 117 or dic.label_id.nunique() != 117:
        sys.exit("label dictionary is not the 117-class dictionary")
    nid = dict(zip(dic.structure, dic.label_id.astype(int)))
    G = build_ids(nid)
    P3 = pd.read_csv(OUT_DIR / P3_RUN / "P3_components_20261112_001915.csv")
    VU, VM = P3[P3.status == "V_U"], P3[P3.status == "V_M"]
    if len(VU) != 217 or len(VM) != 2632:
        sys.exit(f"expected 217 V_U and 2,632 V_M components, found {len(VU)} / {len(VM)}")
    P1C = pd.read_csv(OUT_DIR / P1_RUN / "P1_components_20261111_233936.csv")
    P1C = P1C[P1C.connectivity == CONN].set_index(["case_id", "component_id"])
    P1P = pd.read_csv(OUT_DIR / P1_RUN / "P1_patients_20261111_233936.csv"); P1P = P1P[P1P.connectivity == CONN].set_index("case_id")
    P2C = pd.read_csv(OUT_DIR / P2_RUN / "P2_components_20261111_235603.csv").set_index(["case_id", "component_id"])
    PL = pd.read_csv(OUT_DIR / "SAP_R1_plane_table_20261111_214952.csv").set_index("case_id")
    GR = pd.read_csv(OUT_DIR / P2B_RUN / "P2b_psma_pet_grid_20261112_000903.csv").set_index("case_id")
    cases = sorted(set(VU.case_id) | set(VM.case_id))
    if a.preflight:
        cases = sorted(VU.case_id.unique())[:2]
    run = datetime.now().strftime("%Y%m%d_%H%M%S")
    odir = OUT_DIR / (f"AD_PREFLIGHT_{run}" if a.preflight else f"AD_{run}"); odir.mkdir()
    rows, pairs, imgs = [], [], {}
    for n, cid in enumerate(cases, 1):
        r_, p_ = process(cid, VU[VU.case_id == cid], VM[VM.case_id == cid], G, nid, P1C, P1P, P2C, PL, GR, imgs)
        rows += r_; pairs += p_
        print(f"[{n}/{len(cases)}] {cid}: V_U {len(r_)}; matched pairs {len(p_)} ({sum(x['flagged'] for x in p_)} flagged)  "
              f"({time.time() - t0:.0f} s)", flush=True)
    C = pd.DataFrame(rows); PR = pd.DataFrame(pairs)
    vu_set = VU[VU.case_id.isin(cases)]
    if len(C) != len(vu_set) or abs(C.volume_native_mL.sum() - vu_set.volume_native_mL.sum()) > 1e-9:
        sys.exit("reconciliation failed: distribution does not cover the primary V_U set exactly once")
    if not a.preflight and (len(C) != 217 or C.groupby(["case_id", "component_id"]).size().max() != 1):
        sys.exit("reconciliation failed: 217 components, each once")
    # category table
    cat = []
    for (c, s), d in C.groupby(["category", "subcategory"]):
        cat.append(dict(category=c, subcategory=s, **cell(d)))
    for c, d in C.groupby("category"):
        cat.append(dict(category=c, subcategory="(all)", **cell(d)))
    cat.append(dict(category="total", subcategory="", **cell(C)))
    CAT = pd.DataFrame(cat).sort_values(["category", "subcategory"])
    # flag summary
    fs = []
    if len(PR):
        comp = PR.groupby(["case_id", "fdg_component"]).agg(volume_native_mL=("volume_native_mL", "first"),
                                                           any_flag=("flagged", "any"), any_eval=("evaluable", "any")).reset_index()
        fs = [dict(item="matched pairs", n=len(PR)), dict(item="pairs evaluable", n=int(PR.evaluable.sum())),
              dict(item="pairs flagged", n=int(PR.flagged.sum())),
              dict(item="matched FDG components", n=len(comp), native_mL=float(comp.volume_native_mL.sum()), patients=comp.case_id.nunique()),
              dict(item="FDG components with >= 1 flagged pair", n=int(comp.any_flag.sum()),
                   native_mL=float(comp[comp.any_flag].volume_native_mL.sum()), patients=comp[comp.any_flag].case_id.nunique()),
              dict(item="FDG components with no evaluable pair", n=int((~comp.any_eval).sum()),
                   native_mL=float(comp[~comp.any_eval].volume_native_mL.sum()), patients=comp[~comp.any_eval].case_id.nunique())]
        for (x, y), d in PR[PR.flagged].groupby(["fdg_compartment", "psma_compartment"]):
            fs.append(dict(item=f"flagged pairs: FDG {x} / PSMA {y}", n=len(d)))
    FS = pd.DataFrame(fs)
    out = {}
    for nm, df in (("components", C), ("categories", CAT), ("flag_pairs", PR), ("flag_summary", FS)):
        p = odir / f"AD_{nm}_{run}.csv"; df.to_csv(p, index=False); out[p.name] = md5(p)
    S = [f"ANATOMICAL DISTRIBUTION (SAP 8.3) of the primary unmatched FDG burden: {len(C)} components, "
         f"{C.volume_native_mL.sum():.1f} mL, {C.case_id.nunique()} patients", "category | subcategory | components | native mL | patients"]
    S += [f"  {r.category} | {r.subcategory} | {r.components} | {r.native_mL:.1f} | {r.patients}" for r in CAT.itertuples()]
    S += ["Soft tissue is described by location relative to predefined nodal regions; it is not classified as nodal disease.",
          "Supradiaphragmatic = above the locked T12/L1 plane (operational proxy, not the anatomical diaphragm).", "",
          "ANATOMICAL CONSISTENCY FLAG (separate descriptive audit of matched pairs; not part of the unmatched distribution):"]
    S += [f"  {r.item}: {r.n}" + (f"; {r.native_mL:.1f} mL; {int(r.patients)} patients" if "native_mL" in FS and pd.notna(r.native_mL) else "")
          for r in FS.itertuples()] if len(FS) else ["  no matched pairs in this run"]
    p = odir / f"AD_summary_{run}.txt"; p.write_text("\n".join(S) + "\n", encoding="utf-8"); out[p.name] = md5(p)
    import SimpleITK as sitk, scipy
    man = dict(stage="AD anatomical distribution (SAP 8.3) and consistency flag", mode="PREFLIGHT" if a.preflight else "FULL",
               run=run, script=Path(__file__).name, script_md5=md5(Path(__file__)),
               inputs={f"{d / n}": m for (d, n), m in INPUTS.items()}, image_md5=imgs, cases=cases,
               software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                             scipy=scipy.__version__, SimpleITK=sitk.__version__),
               outputs=out, runtime_s=round(time.time() - t0), status="COMPLETE")
    (odir / f"AD_manifest_{run}.json").write_text(json.dumps(man, indent=2, default=str), encoding="utf-8")
    print("\n".join(S)); print(f"\nOutput: {odir}")


if __name__ == "__main__":
    main()

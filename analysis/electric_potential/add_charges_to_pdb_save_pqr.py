#!/usr/bin/env python3

#################################################################################
#written by GPT-5 with the following prompt:
# Write a python script to take a gromacs topol.top file which references 
# several .itp files via include statements, and a .pdb file describing the same system, 
# and set the charges in the .pdb file using the information from the .top and .itp files.

#and the follow up prompt: 
# yes please match atoms by name
#################################################################################

"""
set_pdb_charges.py

Read a GROMACS .top file (and its included .itp files), extract atom charges
from [ atoms ] sections, and write those charges into a PDB file by matching
on (residue name, atom name). Charges can be written to the occupancy or
B-factor column.

Assumptions:
- GROMACS [ atoms ] lines follow the usual layout:
  nr type resnr residu atom cgnr charge mass [ ... ]
- Residue names in the topology match PDB residue names (e.g., HIS/HID/HIE/HIP,
  NTER/CTER variants, ligand 3-letter codes). If two different molecule types
  use the same residu but with different charges, the first seen wins (a warning
  is printed). You can adjust precedence if needed.
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict

INCLUDE_RE = re.compile(r'^\s*#\s*include\s+"([^"]+)"')
SECTION_RE = re.compile(r'^\s*\[\s*(\w+)\s*\]')
COMMENT_RE = re.compile(r'\s*;.*$')

def _strip_comment(line: str) -> str:
    return COMMENT_RE.sub('', line).rstrip('\n')

def parse_atoms_from_text(text: str):
    """
    Parse all [ atoms ] sections in a given GROMACS topology text.
    Returns a list of dicts with keys: resname, atomname, charge.
    """
    atoms = []
    in_atoms = False
    for raw in text.splitlines():
        line = _strip_comment(raw).strip()
        if not line:
            continue
        m = SECTION_RE.match(line)
        if m:
            in_atoms = (m.group(1).lower() == 'atoms')
            continue
        if not in_atoms:
            continue
        # within [ atoms ]
        # Skip pure comment/empty handled above
        parts = line.split()
        # Expected minimal length: 7 (index 0..6), often 8 with mass
        # Layout: nr type resnr residu atom cgnr charge mass
        if len(parts) < 7:
            continue
        try:
            #resnumber = parts[2]
            resname = parts[3]
            atomname = parts[4]
            charge = float(parts[6])
        except Exception:
            continue
        #JHB WARNING: the index 0:4 is a kluge to make TLCL2 map to TLCL
        atoms.append({'resname': resname[0:3], 'atomname': atomname, 'charge': charge}) #, 'resi': resnumber
    return atoms


def read_topology_recursive(top_path: Path):
    """
    Read a .top or .itp file and all recursively included files.
    Returns concatenated text and a set of visited paths to avoid cycles.
    """
    visited = set()
    texts = []

    def _read(p: Path):
        rp = p.resolve()
        if rp in visited:
            return
        visited.add(rp)
        try:
            content = rp.read_text()
        except Exception as e:
            print(f"Warning: cannot read {rp}: {e}")
            return
        texts.append(content)

        # find includes relative to current file
        for line in content.splitlines():
            m = INCLUDE_RE.match(line)
            if m:
                inc = m.group(1)
                inc_path = (rp.parent / inc).resolve()
                if inc_path.exists():
                    _read(inc_path)
                else:
                    print(f"Warning: included file not found: {inc_path}")

    _read(top_path)
    return "\n".join(texts), visited


def build_residue_templates(top_text: str):
    """
    From concatenated topology text, construct per-residue templates:
    {resname -> {atomname -> charge}}
    If the same (resname, atomname) appears with different charges across files,
    the first occurrence is kept and a warning is printed for conflicts.
    """
    records = parse_atoms_from_text(top_text)
    #print(records)

    templates = defaultdict(dict)
    seen = set()
    for r in records:
        key = (r['resname'], r['atomname'])
        if key in seen:
            # possible duplicate; warn only if charge differs
            old = templates[r['resname']][r['atomname']]
            if abs(old - r['charge']) > 1e-6:
                print(r)
                print(f"Warning: conflicting charge for {key}: keeping {old}, ignoring {r['charge']}")

            continue
        templates[r['resname']][r['atomname']] = r['charge']
        seen.add(key)
    return templates


def parse_pdb_atom_fields(line: str):
    return {
        'atomname': line[12:16].strip(),
        'resname': line[17:20].strip(),
        'chain': line[21:22],
        'resid': line[22:26].strip(),
        'x': float(line[30:38]),
        'y': float(line[38:46]),
        'z': float(line[46:54]),
        'record': line[0:6].strip(),
        'serial': int(line[6:11]),
        'element': line[76:78].strip() if len(line) >= 78 else '',
    }


def remove_terminal_resname(residue):
    if residue == "VAc":
        return "VAL"
    elif residue == "SEn":
        return "SER"
    else:
        return residue


def write_pqr(in_pdb: Path, out_pqr: Path, res_templates: dict, default_radius: float):
    updated = 0
    missing_res = set()
    missing_atom = set()
    with in_pdb.open() as fin, out_pqr.open('w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 54:
                info = parse_pdb_atom_fields(line)
                resmap = res_templates.get(info['resname'])
                if resmap is None:
                    missing_res.add(info['resname'])
                    continue
                charge = resmap.get(info['atomname'])
                if charge is None:
                    missing_atom.add((info['resname'], info['atomname']))
                    continue

                fout.write(
                    f"{info['record']:<6}{info['serial']:>5}  "+
                    f"{info['atomname']:<4}{remove_terminal_resname(info['resname']):>3} "+
                    f"{info['chain']}{info['resid']:>4}    "+
                    f"{info['x']:8.3f}{info['y']:8.3f}{info['z']:8.3f}"+
                    f" {charge:8.4f}{default_radius:8.4f}\n"
                )
                updated += 1
            else:
                # skip non-ATOM/HETATM records in PQR
                continue
    print(f"Updated {updated} atoms in PQR.")
    if missing_res:
        print(f"Residues not in topology: {sorted(missing_res)}")
    if missing_atom:
        print(f"Missing atom matches: {sorted(missing_atom)[:10]}{' ...' if len(missing_atom) > 10 else ''}")


# def parse_pdb_atom_fields(line: str):
#     """
#     Parse minimal fields from a PDB ATOM/HETATM line.
#     Returns dict with indices and values needed; relies on fixed columns.
#     """
#     # PDB format columns (1-based):
#     #  1-6  Record name "ATOM  " / "HETATM"
#     # 13-16 Atom name
#     # 17    altLoc
#     # 18-20 resName
#     # 22    chainID
#     # 23-26 resSeq
#     # 27    iCode
#     # 31-38 x, 39-46 y, 47-54 z
#     # 55-60 occupancy
#     # 61-66 tempFactor (B-factor)
#     # 77-78 element, 79-80 charge (formal)
#     name  = line[12:16].strip()
#     resn  = line[17:21].strip()
#     chain = line[21:22]
#     resi  = line[22:26].strip()
#     return {
#         'atomname': name,
#         'resname': resn,
#         'chain': chain,
#         'resid': resi,
#     }


# def write_pdb_with_charges(in_pdb: Path, out_pdb: Path, res_templates: dict, field: str):
#     """
#     Write charges into the chosen field: 'occupancy' or 'bfactor'.
#     - occupancy columns: 55-60 (0-based slice [54:60])
#     - bfactor columns:   61-66 (0-based slice [60:66])
#     Leaves the other field untouched.
#     """
#     if field not in ('occupancy', 'bfactor'):
#         raise ValueError("field must be 'occupancy' or 'bfactor'")

#     updated = 0
#     missing_res = set()
#     missing_atom = set()

#     with in_pdb.open() as fin, out_pdb.open('w') as fout:
#         for line in fin:
#             if line.startswith(('ATOM', 'HETATM')) and len(line) >= 66:
#                 info = parse_pdb_atom_fields(line)
#                 resmap = res_templates.get(info['resname'])
#                 if resmap is None:
#                     missing_res.add(info['resname'])
#                     fout.write(line)
#                     continue
#                 charge = resmap.get(info['atomname'])
#                 if charge is None:
#                     # Some topologies use names like H1/H2/H3 vs H
#                     # You can add normalization rules here if needed
#                     missing_atom.add((info['resname'], info['atomname']))
#                     fout.write(line)
#                     continue

#                 # Format charge into chosen column
#                 if field == 'occupancy':
#                     # columns 55-60
#                     new = line[:54] + f"{charge:6.3f}" + line[60:]
#                 else:  # bfactor
#                     # columns 61-66
#                     new = line[:60] + f"{charge:6.3f}" + line[66:]
#                 fout.write(new)
#                 updated += 1
#             else:
#                 fout.write(line)

#     print(f"Updated {updated} atoms with charges into {field}.")
#     if missing_res:
#         print(f"Note: {len(missing_res)} PDB residue names not found in topology: {sorted(missing_res)}")
#     if missing_atom:
#         samples = sorted(list(missing_atom))[:10]
#         extra = f" (+{len(missing_atom)-10} more)" if len(missing_atom) > 10 else ""
#         print("Note: Missing atom names for residue/atom pairs (showing up to 10): "
#               f"{samples}{extra}")


def main():
    ap = argparse.ArgumentParser(description="Set PDB charges from GROMACS .top/.itp files by matching (resname, atomname).")
    ap.add_argument('-t', '--top', required=True, help='Path to topol.top')
    ap.add_argument('-p', '--pdb', required=True, help='Input PDB file')
    ap.add_argument('-o', '--out', required=True, help='Output PDB file')
    ap.add_argument('--radius', type=float, default=1.5, help='Default atomic radius (Å) if none provided')
    # ap.add_argument('--field', choices=['occupancy', 'bfactor'], default='bfactor',
    #                 help='Which PDB field to write charges into (default: bfactor)')
    args = ap.parse_args()

    top_path = Path(args.top)
    pdb_in = Path(args.pdb)
    pdb_out = Path(args.out)

    top_text, visited = read_topology_recursive(top_path)
    print(f"Parsed topology files: {len(visited)}")
    templates = build_residue_templates(top_text)
    print(f"Built residue templates for {len(templates)} residue names.")
    write_pqr(pdb_in, pdb_out, templates, args.radius)
    print(f"Wrote {pdb_out}")


if __name__ == '__main__':
    main()
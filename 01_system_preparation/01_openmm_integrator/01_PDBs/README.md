### Author: William Glass

## Notes on preparing the ACE2:RBD Complex using ISOLDE

The fully glycosylated ACE2:RBD complex was constructed by aligning the previously constructed ACE2 and RBD systems with the refined `6m0j` [structure](https://github.com/thorn-lab/coronavirus_structural_task_force/tree/master/pdb/surface_glycoprotein/SARS-CoV-2/6m0j).

It was immediately clear that there were sever clashes between residues at the interface of ACE2 and the RBD. In order to accurately describe the interface [ISOLDE](https://isolde.cimr.cam.ac.uk/) was used for structural refinement.

### Steps taken in ISOLDE

Below describes the steps taken to refine the ACE2:RBD complex.

1. After evaluation of the initial aligned ACE2:RBD complex (`0_RBD_ACE2_complex_from_1r42ref_6m0j_FyllGlycos_alignto6m0jRefined_NOCAPS.pdb`) in ChimeraX and ISOLDE, it was clear there were difficulties in assigning the correct hydrogens to glycans (mainly due to the manual fitting of glycans to each glycosylation site). This was solved by [Tristan Croll](https://github.com/tristanic), where a "hydrogens corrected" system was generated: `1_refining_ace2_rbd_complex_hydrogens_corrected.pdb`.

2. Within the "hydrogens corrected" model, glycans had been renamed to PDB format (e.g. NAG, MAN, etc...) and bonding between glycans had already been specified. These changes were applied using the scripts detailed at the end of this README.

3. In ChimeraX, the refined `6m0j` CST structure was loaded in alongside the `1_refining_ace2_rbd_complex_hydrogens_corrected.pdb` model. These were then aligned:
    * Command line: `match #2 to #1`

4. All crystal waters were removed due to clashes:
    * Command line: `delete :HOH`

5. Launch ISOLDE:
    * Command line: `isolde start`

6. Within ISOLDE it is possible to view clashes between residues (ISOLDE -> Validate- > Clashes). The resulting list was used to then correct any major clashes (minor clashes / very close residues can be handled). Correction was performed via selection of the alpha carbon (`ctrl` + `click`) and using the rotamer selection tool. Rotamers were chosen to best match those present in the refined `6m0j` structure (that is shown whilst using ISOLDE).
    * Rotamers changed:
        * MET82 (ACE2) -> it was clashing with F486 (RBD)
        * K31 (ACE2) - > it was clashing with T489 (RBD)
        * Flipped V107 (ACE2) as to match 6m0j better and reduce clash with the nearby N103 glycan.

7. Two sets of [custom scripts](https://github.com/tristanic/isolde/blob/75f3ad7176a8fc200681ceb1b9197562f732ba5e/isolde/src/atomic/building/build_utils.py#L259) were used to create disulphide bonds within both ACE2 and RBD chains:
    ```python
    current, possible, ambiguous = current_and_possible_disulfides(m)
    for cys_pair in possible:
        create_disulfide(*cys_pair)
    ```
Each bond was checked and if any had not been bonded a custom `make_bonds.py` script was used to enforce bonding and any remaining hydrogens removed.

8. After assigning new rotameric states to residues the torsion restraints need to be reapplied. This was achieved using:
    * Command line: `isolde restrain torsions #1/A templ #2/A angleRange 180 springConstant 1000` and `isolde restrain torsions #1/E templ #2/E angleRange 180 springConstant 1000`

9. ISOLDE simulations were run at interface residues until they match roughly the refined `6m0j` structure residue rotamers.

10. Glycans were then renamed to their original format (see script below). This also renamed ASN to NLN and CYS to CYX where appropriate for the AMBER `tleap` input.

11. The refined `6m0j` structure was removed and the final refined ACE2:RBD saved as `2_refined_ace2_rbd_complex.pdb`

12. Hydrogens were removed from the final structure due to naming issues (where `tleap` added extra hydrogens). By removing these, `tleap` subsequently added hydrogens correctly. Secondary structure information was removed (as this is created in `tleap` later). Additionally, `CONECT` records were removed as these are not read by `tleap` and when adding termini caps later, will not correspond to the correct residues.

* **Note**: investigating how `tleap` can be used with `CONECT` records would reduce bottlenecks in the described pipeline above.

---

Scripts used within ChimeraX / ISOLDE

* Renaming glycans to PDB format (NAG, MAN, etc)
    ```python
    def glycam_to_pdb_residue_names(model):

        from chimerax.atomic import Residue

        non_protein = model.residues[model.residues.polymer_types!=Residue.PT_AMINO]

        from chimerax.isolde.openmm.amberff import glycam

        def is_glycam(rname):

            return (rname[1:] in glycam.glycam_suffix_to_ccd_name.keys() and
                rname[0] in glycam._glycam_prefix.values())

        glycam_name_map = {r: r.name for r in non_protein if is_glycam(r.name)}
        
        import numpy

        anchor_names = list(glycam._anchor_name_map.values())

        
        anchor_names = anchor_names + ['N'+aname for aname in anchor_names] + ['C'+aname for aname in anchor_names]
        
        anchors = m.residues[numpy.in1d(m.residues.names, anchor_names)]
        
        anchor_name_map = {r: r.name for r in anchors}
        
        for r, gname in glycam_name_map.items():
            r.name = glycam.glycam_suffix_to_ccd_name[r.name[1:]]
        glycam_to_pdb_anchor = {val: key for key,val in glycam._anchor_name_map.items()}
        
        for r, gname in anchor_name_map.items():
            r.name = glycam_to_pdb_anchor[r.name[-3:]]
        glycam_name_map.update(anchor_name_map)
        
        return glycam_name_map

    def revert_to_glycam_names(glycam_name_map):
        
        for r, rname in glycam_name_map.items():
            r.name = rname
    ```

* Converting glycan names back from PDB format
    ```python
    from chimerax.isolde.openmm.openmm_interface import find_residue_templates
    ff = session.isolde.forcefield_mgr['amber14']
    template_dict = find_residue_templates(m.residues, ff)
    for i, template_name in template_dict.items():
        if 'GLYCAM' in template_name:
            m.residues[i].name = template_name.split('_')[1]
        elif '_' not in template_name:
            m.residues[i].name = template_name
    ```
* Make bonds script (located in a startup folder for general usage)
    ```python
    # @Author: Tristan Croll <tic20>
    # @Date:   07-Apr-2020
    # @Email:  tic20@cam.ac.uk
    # @Last modified by:   tic20
    # @Last modified time: 07-Apr-2020
    # @License: Lesser GNU Public License version 3.0 (see LICENSE.md)
    # @Copyright: 2016-2019 Tristan Croll
    from chimerax.atomic import selected_atoms
    from chimerax.core.errors import UserError
    sel = selected_atoms(session)
    if len(sel) != 2:
        raise UserError('Must have exactly two atoms selected!')
    us = sel.unique_structures
    if len(us) != 1:
        raise UserError('Both atoms must be from the same structure!')
    m = us[0]
    m.new_bond(*sel)
    ```

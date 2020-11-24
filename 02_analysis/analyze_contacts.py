import argparse
from multiprocessing import Manager, Pool
from pathlib import Path

import MDAnalysis as mda
import MDAnalysis.transformations as trans
import numpy as np
import pandas as pd
from MDAnalysis.analysis import align, distances, rms
from MDAnalysis.analysis.hbonds import WaterBridgeAnalysis

parser = argparse.ArgumentParser(description="Process FAH data.")
parser.add_argument(
    "-pdb_file",
    dest="pdb_file",
    type=str,
    help="the equilibrated structure, with canonical residue numbering",
)
parser.add_argument(
    "-traj_path",
    dest="traj_path",
    type=str,
    help="the location of the FAH trajectories i.e. /path/to/traj/PROJXXXXX/RUNX",
)
parser.add_argument(
    "-stride",
    dest="stride",
    type=int,
    help="the stride length. This will be the number of CLONES skipped during analysis",
    default=None,
)
parser.add_argument(
    "-fah_project_code",
    dest="project_code",
    type=str,
    help="the FAH project code, numbers only",
)
parser.add_argument(
    "-split",
    dest="split_value",
    type=int,
    help="the number used to split the list of trajectories into smaller chunks",
    default=100,
)
parser.add_argument(
    "-n_clones",
    dest="n_clones",
    type=int,
    help="the total number of clones to analyse",
    default=10000,
)
parser.add_argument(
    "-n_gens",
    dest="n_gens",
    type=int,
    help="the total number of gens (per CLONE) to analyse",
    default=1000,
)
parser.add_argument(
    "-mutant",
    dest="mutant",
    type=str,
    choices=["WT"],
    help="the mutant to analyse e.g. WT",
)
parser.add_argument(
    "-starting_gen",
    dest="starting_gen",
    type=int,
    help="GEN number to start analysis from i.e. if 5 chosen, only use trajectories >= GEN5 in each CLONE",
    default=0,
)
# TODO define -clone_file type
parser.add_argument(
    "-clone_file",
    dest="clone_file",
    help="a .txt file containing the CLONES to analyse on each line. Only these CLONE numbers will be analysed if parsed.",
    default=None,
)
parser.add_argument(
    "-stride_frames",
    dest="stride_frames",
    type=int,
    help="the number of frames to stride in a GEN",
    default=1,
)

args = parser.parse_args()

# define the class to help group the two chains together
# The GroupHug class was created by Richard Gowers (https://github.com/richardjgowers)
# in response to this question on the MDAnalysis forum:
# https://groups.google.com/forum/#!topic/mdnalysis-discussion/umDpvbCmQiE
class GroupHug:
    def __init__(self, center, *others):
        self.c = center
        self.o = others

    @staticmethod
    def calc_restoring_vec(ag1, ag2):
        box = ag1.dimensions[:3]
        dist = ag1.center_of_mass() - ag2.center_of_mass()

        return box * np.rint(dist / box)

    def __call__(self, ts):
        # loop over other atomgroups and shunt them into nearest image to center
        for i in self.o:
            rvec = self.calc_restoring_vec(self.c, i)

            i.translate(+rvec)

        return ts


# define the function to collect RMSDs
def populate_dict(chunk, dictionary, pdb_file, mutant_sel, project_code, frames_to_stride):

    interface_selection_strings = {
        "rbd": "segid A and (backbone and (resid 403 or resid 417 or resid 439 or resid 445-447 or resid 449 or resid 453 or resid 455 or resid 456 or resid 473-477 or resid 484-487 or resid 489 or resid 490 or resid 493-503 or resid 505 or resid 506))",
        "ace2": "segid C and (backbone and (resid 18 or resid 21 or resid 23-32 or resid 33-39 or resid 41 or resid 42 or resid 45 or resid 75 or resid 76 or resid 78-84))",
        "rbd_and_ace2": "(segid A and (backbone and (resid 403 or resid 417 or resid 439 or resid 445-447 or resid 449 or resid 453 or resid 455 or resid 456 or resid 473-477 or resid 484-487 or resid 489 or resid 490 or resid 493-503 or resid 505 or resid 506))) or (segid C and (backbone and (resid 18 or resid 21 or resid 23-32 or resid 33-39 or resid 41 or resid 42 or resid 45 or resid 75 or resid 76 or resid 78-84)))",
    }

    # Create a dictionary containing selection strings for MDAnalysis
    # residues 417 and 439 are not named since these are mutated across systems
    # segid C = ACE2, segid A = RBD

    # TODO remove project keys, not sure they are needed
    proj_mutant_dict = {
        "17311": {
            "WT": {
                "D30": "segid C and (resid 30 and name OD1 OD2)",
                "res417": "segid A and (resid 417 and name NZ)",
                "E329": "segid C and (resid 329 and name OE1 OE2)",
                "res439": "segid A and (resid 439 and name ND2)",
                "K31": "segid C and (resid 31 and name NZ)",
                "E484": "segid A and (resid 484 and name OE1 OE2)",
                "E35": "segid C and (resid 35 and name OE1 OE2)",
                "K31": "segid C and (resid 31 and name NZ)",
                "Q493": "segid A and (resid 493 and name NE2 OE1)",
                "K353": "segid C and (resid 353 and name NZ)",
                "G496bb": "segid A and (resid 496 and name O C CA N)",
                "D38": "segid C and (resid 38 and name OD1 OD2)",
                "Y449": "segid A and (resid 449 and name CG CD1 CE1 CZ CE2 CD2 OH)",
                "Q42": "segid C and (resid 42 and name NE2 OE1)",
                "K353bb": "segid C and (resid 353 and name O C CA N)",
                "G502bb": "segid A and (resid 502 and name O C CA N)",
            },
        },
    }

    # set the reference to be the equilibrated structure
    ref = mda.Universe(pdb_file)
    ref.trajectory[0]  # there is only one frame anyway, but just to be sure
    ref_bb = ref.select_atoms("backbone")  # the ref for RMSD calcs later

    # set reference interfaces
    ref_rbd_interface_bb = ref.select_atoms(interface_selection_strings["rbd"])

    ref_ace2_interface_bb = ref.select_atoms(interface_selection_strings["ace2"])

    ref_whole_interface_bb = ref.select_atoms(
        interface_selection_strings["rbd_and_ace2"]
    )

    reference = ref.select_atoms(
        "not resname Na+ Cl- HOH"
    )  # the ref for transforms later

    for traj in chunk:
        print("--> Analysing trajectory: ", traj)

        mobile = mda.Universe(pdb_file, traj)
        # centre the two protein chains in the box
        # this stops chains jumping across PBC
        chainA = mobile.select_atoms("segid A or segid B")  # RBD + glycans
        chainB = mobile.select_atoms("segid C or segid D")  # ACE2 + glycans
        ions = mobile.select_atoms("resname Na+ Cl- HOH")
        protein = mobile.select_atoms("not resname Na+ Cl- HOH")

        transforms = [
            trans.unwrap(protein),
            GroupHug(chainA, chainB),
            trans.center_in_box(protein, wrap=False, center="geometry"),
            trans.wrap(ions),
            trans.fit_rot_trans(protein, reference),
        ]

        print("--> Centring protein chains in the box")
        mobile.trajectory.add_transformations(*transforms)

        # loop over each frame in the current trajectory, with a defined stride
        for ts in mobile.trajectory[::frames_to_stride]:

            print(ts.frame)

            # calculate the key interactions
            # RBD --- ACE2

            # D30 --- K417
            D30 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["D30"])
            res417 = mobile.select_atoms(
                proj_mutant_dict[project_code][mutant]["res417"]
            )

            D30_res417_dist_mindist = np.min(
                distances.distance_array(D30.positions, res417.positions)
            )

            # E329 --- N439
            E329 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["E329"])
            res439 = mobile.select_atoms(
                proj_mutant_dict[project_code][mutant]["res439"]
            )

            E329_res439_dist_mindist = np.min(
                distances.distance_array(E329.positions, res439.positions)
            )

            # E484 --- K31
            K31 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["K31"])
            E484 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["E484"])

            E484_K31_dist_mindist = np.min(
                distances.distance_array(E484.positions, K31.positions)
            )

            # E35 --- K31
            E35 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["E35"])

            E35_K31_dist_mindist = np.min(
                distances.distance_array(E35.positions, K31.positions)
            )

            # E35 --- Q493
            Q493 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["Q493"])

            E35_Q493_dist_mindist = np.min(
                distances.distance_array(E35.positions, Q493.positions)
            )

            # Additional interactions
            K31_Q493_dist_mindist = np.min(
                distances.distance_array(K31.positions, Q493.positions)
            )

            # K353 --- G496 (K353 to G496 backbone)
            K353 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["K353"])
            G496bb = mobile.select_atoms(
                proj_mutant_dict[project_code][mutant]["G496bb"]
            )

            K353_G496bb_dist_mindist = np.min(
                distances.distance_array(K353.positions, G496bb.positions)
            )

            # D38 --- Y449
            D38 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["D38"])
            Y449 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["Y449"])

            D38_Y449_dist_mindist = np.min(
                distances.distance_array(D38.positions, Y449.positions)
            )

            # Q42 --- Y449
            Q42 = mobile.select_atoms(proj_mutant_dict[project_code][mutant]["Q42"])

            Q42_Y449_dist_mindist = np.min(
                distances.distance_array(Q42.positions, Y449.positions)
            )

            # K353bb --- G502bb (Backbone to backbone)
            K353bb = mobile.select_atoms(
                proj_mutant_dict[project_code][mutant]["K353bb"]
            )
            G502bb = mobile.select_atoms(
                proj_mutant_dict[project_code][mutant]["G502bb"]
            )

            K353bb_G502bb_dist_mindist = np.min(
                distances.distance_array(K353bb.positions, G502bb.positions)
            )

            # sort out names for the dict
            traj_split = traj.split("/")
            key_name = f"{traj_split[6]}/{traj_split[7]}/{traj_split[8]}_{ts.frame}"

            # populate the dict - with placeholder keys for now
            dictionary[key_name] = {
                "d30_res417_mindist": D30_res417_dist_mindist,
                "e329_res439_mindist": E329_res439_dist_mindist,
                "e484_k31_mindist": E484_K31_dist_mindist,
                "e35_k31_mindist": E35_K31_dist_mindist,
                "e35_q493_mindist": E35_Q493_dist_mindist,
                "q493_k31_mindist": K31_Q493_dist_mindist,
                "k353_g496bb_mindist": K353_G496bb_dist_mindist,
                "d38_y449_dist_mindist": D38_Y449_dist_mindist,
                "q42_y449_dist_mindist": Q42_Y449_dist_mindist,
                "k353bb_g502bb_dist_mindist": K353bb_G502bb_dist_mindist,
            }


if __name__ == "__main__":

    # Define everything
    pdb_file = args.pdb_file
    mutant = args.mutant
    project_code = args.project_code
    traj_path = args.traj_path
    n_clones = args.n_clones  # Total No. CLONES = 10000
    n_gens = args.n_gens
    stride = args.stride
    clone_file = args.clone_file
    starting_gen = args.starting_gen
    stride_frames = args.stride_frames

    # Print for a sanity check
    print(f"--> Using PDB file: {pdb_file}")
    print(f"--> Analysing {mutant} system")
    print(f"--> Using project code: {project_code}")

    # read in the trajs as a list
    print("\n--> Creating list of trajectories")
    print(f"--> Starting GEN number: {starting_gen}")
    print(f"--> Number of GENS to check: {n_gens}")
    print(f"--> Number of frames to stride per GEN: {stride_frames}")

    traj_list = []

    # TODO clean up all of this
    if clone_file is not None: # if we are using a input file of clones
        with open(clone_file) as f:
            temp_file = f.readlines()
            clones_to_analyse = [line.rstrip("\n") for line in temp_file]

        if stride is not None:
            for clone_number in clones_to_analyse[::stride]:
                for j in range(starting_gen, n_gens):

                    xtc_file = (
                        f"{traj_path}/CLONE{clone_number}/results{j}/positions.xtc"
                    )
                    xtc_file_path = Path(xtc_file)

                    if xtc_file_path.is_file():
                        traj_list.append(xtc_file)
        else:
            for clone_number in clones_to_analyse:
                for j in range(starting_gen, n_gens):

                    xtc_file = (
                        f"{traj_path}/CLONE{clone_number}/results{j}/positions.xtc"
                    )
                    xtc_file_path = Path(xtc_file)

                    if xtc_file_path.is_file():
                        traj_list.append(xtc_file)

    else: # if we are running without an input file
        for i in range(n_clones):

            # TODO clean this up, it's horrible
            if stride is not None:  # stride over CLONES
                if i % stride == 0:
                    for j in range(starting_gen, n_gens):

                        xtc_file = f"{traj_path}/CLONE{i}/results{j}/positions.xtc"
                        xtc_file_path = Path(xtc_file)

                        if xtc_file_path.is_file():
                            traj_list.append(xtc_file)
            else:  # don't stride, read all clones
                for j in range(starting_gen, n_gens):

                    xtc_file = f"{traj_path}/CLONE{i}/results{j}/positions.xtc"
                    xtc_file_path = Path(xtc_file)

                    if xtc_file_path.is_file():
                        traj_list.append(xtc_file)

    # break the traj list up ready for multiprocessing
    print("--> Breaking trajectory list into smaller chunks...")
    split_value = args.split_value
    traj_chunks = [
        traj_list[x : x + split_value] for x in range(0, len(traj_list), split_value)
    ]

    print("\n--> Staring pool...")
    manager = Manager()
    cluster_dict = manager.dict()  # allows the dict to be accessed across processes
    pool = (
        Pool()
    )  # use all available cores, otherwise specify the number you want as an argument
    for chunk in traj_chunks:
        pool.apply_async(
            populate_dict,
            args=(chunk, cluster_dict, pdb_file, mutant, project_code, stride_frames),
        )
    pool.close()
    pool.join()

    # create the dictionary / CSV / DataFrame to analyse later
    print("--> Multiprocessing complete, creating CSV")
    d = cluster_dict.copy()

    df = pd.DataFrame(d)
    df = df.T
    df.reset_index(inplace=True)

    if stride is not None:
        df.to_csv(
            f"PROJ{project_code}_gathered_strided{stride}clones_rbd_ace2_data_RMSD_SALTBRIDGES_{mutant}.csv"
        )
    else:
        df.to_csv(
            f"PROJ{project_code}_gathered_rbd_ace2_data_RMSD_SALTBRIDGES_{mutant}.csv"
        )

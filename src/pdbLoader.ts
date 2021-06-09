import * as pdb from "./libs/pdb";
import * as base from "./libs/base";
import * as THREE from "three";

function loadPDB(source_file: string, strand_dir, models_as_strands) {
    let pdb_strands = [];
    let pdb_strand = [];
    let old_residue: number;
    let old_chain: string;

    let nn: pdb.Nucleotide;
    for (let line of source_file.split("\n")) {
        line = line.trim();
        let na = new pdb.Atom(line);
        if (line.startsWith("ATOM")) {
            if (old_chain !== "") {
                if (na.chain_id !== old_chain && pdb_strand.length !== 0) {
                    base.Logger.log(
                        `WARNING: a TER statement separating different strands (${
                            na.chain_id
                        } and ${old_chain}) is missing`
                    );
                    pdb_strands.push(pdb_strand);
                    pdb_strand = [];
                } else if (
                    na.chain_id === old_chain &&
                    pdb_strand.length === 0
                ) {
                    base.Logger.log(
                        `WARNING: a TER statement separates strands having the same chain id (${
                            na.chain_id
                        })`
                    );
                }
            }
            if (na.alternate !== "") {
                if (na.alternate === "A" || na.alternate === "1") {
                    base.Logger.log(
                        `Alternate location for atom '${na.name}' of residue '${
                            na.residue
                        }' encountered, using the line marked with the '${
                            na.alternate
                        }' character.`
                    );
                }
            }
            if (na.residue_idx !== old_residue) {
                nn = new pdb.Nucleotide(na.residue, na.residue_idx);
                if (strand_dir) {
                    pdb_strand.push(nn);
                } else {
                    pdb_strand.splice(0, 0, nn);
                }
                old_residue = na.residue_idx;
            }
            nn.add_atom(na);
            old_chain = na.chain_id;
        } else if (line.startsWith("MODEL")) {
            if (!models_as_strands) {
                let N_model = line.split(" ")[1];
                base.Logger.log(
                    `MODEL line detected: using the first MODEL encountered (${N_model})`
                );
            }
        } else if (line.startsWith("ENDMDL")) {
            if (!models_as_strands) {
                // by default we treat ENDMDL as the end of the file
                break;
            } else {
                // otherwise we treat it as the end of the strand
                if (pdb_strand.length > 0) {
                    pdb_strands.push(pdb_strand);
                    pdb_strand = [];
                }
            }
        } else if (line.startsWith("TER")) {
            pdb_strands.push(pdb_strand);
            pdb_strand = [];
            // if the file does not contain any TER line we need to manually add the current strand to the list of strands
        } else if (
            line === "END" &&
            pdb_strands.length === 0 &&
            pdb_strand.length > 0
        ) {
            pdb_strands.push(pdb_strand);
            pdb_strand = [];
        }
    }

    // sometimes files just end (without any END or TER line)
    if (pdb_strand.length > 0) {
        pdb_strands.push(pdb_strand);
    }
    let box_low = new THREE.Vector3(1e6, 1e6, 1e6);
    let box_high = new THREE.Vector3(-1e6, -1e6, -1e6);
    for (const nucl of [].concat(...pdb_strands)) {
        let com = nucl.get_com();
        for (let i = 0; i < 3; i++) {
            if (com.getComponent(i) < box_low.getComponent(i)) {
                box_low.setComponent(i, com.getComponent(i));
            } else if (com.getComponent(i) > box_high.getComponent(i)) {
                box_high.setComponent(i, com.getComponent(i));
            }
        }
    }

    let L =
        Math.max(
            ...box_high
                .clone()
                .sub(box_low)
                .toArray()
        ) *
        2 *
        pdb.FROM_ANGSTROM_TO_OXDNA;
    L = Math.ceil(L);
    let box = new THREE.Vector3(L, L, L);

    base.Logger.log(
        `Using a box of size ${L} in oxDNA units (twice as big as the PDB bounding box size)`
    );

    let system = new base.System(box);

    for (const pdb_strand of pdb_strands) {
        let strand = new base.Strand();
        for (const nucl of pdb_strand) {
            nucl.compute_as();
            let com = nucl
                .get_com()
                .clone()
                .multiplyScalar(pdb.FROM_ANGSTROM_TO_OXDNA);
            let new_oxDNA_nucl = new base.Nucleotide(
                com,
                nucl.a1,
                nucl.a3,
                nucl.base[0]
            );
            strand.add_nucleotide(new_oxDNA_nucl);
        }
        system.add_strand(strand);
    }

    return system;
}

export { loadPDB };

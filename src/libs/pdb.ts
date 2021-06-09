import { Logger } from "./base";
import * as THREE from "three";
import * as utils from "./utils";

const BASE_SHIFT = 1.13;
const COM_SHIFT = 0.5;
const FROM_OXDNA_TO_ANGSTROM = 8.518;
const FROM_ANGSTROM_TO_OXDNA = 1 / FROM_OXDNA_TO_ANGSTROM;

const NAME_TO_BASE = new Map([
    ["ADE", "A"],
    ["CYT", "C"],
    ["GUA", "G"],
    ["THY", "T"],
    ["URA", "U"]
]);

const BASES = ["A", "T", "G", "C", "U"];

let RNA_warning_printed = false;
class Nucleotide {
    base: string;
    idx: number;
    name: string;
    base_atoms: Atom[];
    phosphate_atoms: Atom[];
    sugar_atoms: Atom[];
    named_atoms: {};
    ring_names: string[];
    chain_id: string;
    a1: THREE.Vector3;
    a3: THREE.Vector3;
    a2: THREE.Vector3;
    check: number;

    constructor(name: string, idx: number) {
        this.name = name.trim();
        if ([...NAME_TO_BASE.keys()].includes(this.name)) {
            this.base = NAME_TO_BASE[this.name];
        } else if (this.name in BASES) {
            if (this.name == "U" && !RNA_warning_printed) {
                Logger.log(
                    "WARNING: unsupported uracil detected: use at your own risk"
                );
                RNA_warning_printed = true;
            }
            this.base = this.name;
        } else {
            this.base = name.slice(1);
        }
        this.idx = idx;
        this.base_atoms = [];
        this.phosphate_atoms = [];
        this.sugar_atoms = [];
        this.named_atoms = {};
        this.ring_names = ["C2", "C4", "C5", "C6", "N1", "N3"];
    }

    get atoms() {
        return [].concat(
            this.base_atoms,
            this.phosphate_atoms,
            this.sugar_atoms
        );
    }

    add_atom(a) {
        if (a.name.includes("P") || a.name == "HO5'") {
            this.phosphate_atoms.push(a);
        } else if (a.name.includes("'")) {
            this.sugar_atoms.push(a);
        } else {
            this.base_atoms.push(a);
        }
        this.named_atoms[a.name] = a;
        if (this.chain_id === undefined) {
            this.chain_id = a.chain_id;
        }
    }

    get_com(atoms?) {
        if (atoms == undefined) {
            atoms = this.atoms;
        }
        let com = new THREE.Vector3(0, 0, 0);
        for (const a of atoms) {
            com.add(a.pos);
        }
        return com.divideScalar(atoms.length);
    }

    compute_a3() {
        const base_com = this.get_com(this.base_atoms);
        // the O4' oxygen is always (at least for non pathological configurations, as far as I know) oriented 3' -> 5' with respect to the base's centre of mass
        const parallel_to = this.named_atoms["O4'"].pos.clone().sub(base_com);
        this.a3 = new THREE.Vector3(0, 0, 0);

        for (const perm of utils.permutations<string>(this.ring_names, 3)) {
            let p = this.named_atoms[perm[0]];
            let q = this.named_atoms[perm[1]];
            let r = this.named_atoms[perm[2]];
            let v1 = p.pos.clone().sub(q.pos);
            let v2 = p.pos.clone().sub(r.pos);
            v1.divideScalar(Math.sqrt(v1.dot(v1)));
            v2.divideScalar(Math.sqrt(v2.dot(v2)));
            if (Math.abs(v1.dot(v2)) > 0.01) {
                let a3 = v1.clone().cross(v2);
                a3.divideScalar(Math.sqrt(a3.dot(a3)));
                if (a3.dot(parallel_to) < 0) {
                    a3.negate();
                }
                this.a3.add(a3);
            }
        }

        this.a3.divideScalar(Math.sqrt(this.a3.dot(this.a3)));
    }

    compute_a1() {
        let pairs: string[][];
        if (
            this.name.includes("C") ||
            this.name.includes("T") ||
            this.name.includes("U")
        ) {
            pairs = [["N3", "C6"], ["C2", "N1"], ["C4", "C5"]];
        } else {
            pairs = [["N1", "C4"], ["C2", "N3"], ["C6", "C5"]];
        }

        this.a1 = new THREE.Vector3(0, 0, 0);
        for (const pair of pairs) {
            let p = this.named_atoms[pair[0]];
            let q = this.named_atoms[pair[1]];
            let diff = p.pos.clone().sub(q.pos);
            this.a1.add(diff);
        }

        this.a1.divideScalar(Math.sqrt(this.a1.dot(this.a1)));
    }

    compute_as() {
        this.compute_a1();
        this.compute_a3();
        this.a2 = this.a3.clone().cross(this.a1);
        this.check = Math.abs(this.a1.dot(this.a3));
    }

    correct_for_large_boxes(box) {
        // Fix odd js rounding convenions
        let rint = n => Math.sign(n) * Math.round(Math.abs(n));
        for (const atom of this.atoms) {
            atom.shift(-rint(atom.pos / box) * box);
        }
    }

    to_pdb(
        chain_identifier,
        print_H,
        residue_serial,
        residue_suffix,
        residue_type
    ) {
        let res = [];
        let phosphorus: Atom, C: Atom, O3prime: Atom, O5prime: Atom;
        for (const a of this.atoms) {
            if (!print_H && a.name.includes("H")) {
                continue;
            }
            if (residue_type == "5") {
                if ("P" in a.name) {
                    if (a.name == "P") {
                        phosphorus = a;
                    }
                    continue;
                } else if (a.name == "O5'") {
                    C = a;
                }
            } else if (residue_type == "3") {
                if (a.name == "O3'") {
                    O3prime = a;
                }
            }
            res.push(
                a.to_pdb(chain_identifier, residue_serial, residue_suffix)
            );
        }

        // if the residue is a 3' or 5' end, it requires one more hydrogen linked to the O3' or O5', respectively
        if (residue_type == "5") {
            let new_hydrogen = phosphorus.clone();
            new_hydrogen.name = "HO5'";

            // we put the new hydrogen at a distance 1 Angstrom from the O5' oxygen along the direction that, in a regular nucleotide, connects O5' and P
            let dist_P_O = phosphorus.pos.sub(O5prime.pos);
            dist_P_O.multiplyScalar(1 / Math.sqrt(dist_P_O.dot(dist_P_O)));
            new_hydrogen.pos = O5prime.pos.clone().add(dist_P_O);
            res.push(
                new_hydrogen.to_pdb(
                    chain_identifier,
                    residue_serial,
                    residue_suffix
                )
            );
        } else if (residue_type == "3") {
            let new_hydrogen = O3prime.clone();
            new_hydrogen.name = "HO3'";

            // we put the new hydrogen at a distance 1 Angstrom from the O3' oxygen along a direction which is a linear combination of the three
            // orientations that approximately reproduce the crystallographic position
            let new_distance = this.a2
                .clone()
                .multiplyScalar(0.2)
                .sub(this.a1.clone().multiplyScalar(0.2))
                .sub(this.a3);
            new_distance.multiplyScalar(
                1 / Math.sqrt(new_distance.dot(new_distance))
            );
            new_hydrogen.pos = O3prime.pos.clone().add(new_distance);
            res.push(
                new_hydrogen.to_pdb(
                    chain_identifier,
                    residue_serial,
                    residue_suffix
                )
            );
        }

        return res.join("\n");
    }

    to_mgl() {
        let res = [];
        for (const a of this.atoms) {
            res.push(a.to_mgl());
        }

        return res.join("\n");
    }

    rotate(R) {
        let com = this.get_com();
        for (const a of this.atoms) {
            a.pos = R.dot(a.pos.clone().sub(com)).add(com);
        }
        this.compute_as();
    }

    set_com(new_com) {
        let com = this.get_com();
        for (const a of this.atoms) {
            a.pos.add(
                new_com.sub(com).sub(this.a1.clone().multiplyScalar(COM_SHIFT))
            );
        }
    }

    /*
    set_base (new_base_com){
        let atoms = this.ring_names.forEach(k, v in this.named_atoms.items() if k=>v)
        ring_com = this.get_com(atoms)
        for(const a of this.atoms)  {
            a.pos += new_base_com - ring_com - BASE_SHIFT * this.a1
        }

        this.compute_as()
    }
    */
}

class Atom {
    name: string;
    alternate: string;
    residue: string;
    chain_id: string;
    residue_idx: number;
    pos: THREE.Vector3;

    static serial_atom = 1;

    constructor(pdb_line: string) {
        // http://cupnet.net/pdb-format/
        this.name = pdb_line.slice(12, 16).trim();

        // some PDB files have * in place of '
        if (this.name.includes("*")) {
            this.name = this.name.replace("*", "'");
        }

        this.alternate = pdb_line[16];
        this.residue = pdb_line.slice(17, 20).trim();
        this.chain_id = pdb_line.slice(21, 22).trim();
        this.residue_idx = parseInt(pdb_line.slice(22, 26));
        this.pos = new THREE.Vector3(
            parseFloat(pdb_line.slice(30, 38)),
            parseFloat(pdb_line.slice(38, 46)),
            parseFloat(pdb_line.slice(46, 54))
        );
    }

    clone() {
        let c = new Atom("");
        c.name = this.name.slice();
        c.alternate = this.alternate.slice();
        c.residue = this.residue.slice();
        c.chain_id = this.chain_id.slice();
        c.residue_idx = this.residue_idx;
        c.pos = this.pos.clone();
        return c;
    }

    is_hydrogen() {
        this.name.includes("H");
    }

    shift(diff: THREE.Vector3) {
        this.pos.add(diff);
    }

    to_pdb(chain_identifier, residue_serial, residue_suffix) {
        let residue = this.residue + residue_suffix;
        let f = (v, n) => `${v}`.padStart(n, " ");
        let res =
            `ATOM  ${f(Atom.serial_atom, 5)} ${f(this.name, 4)} ${f(
                residue,
                3
            )} ` +
            `${f(chain_identifier, 3)}${f(residue_serial, 1)}    ` +
            `${f(this.pos.x.toFixed(3), 8)}`+
            `${f(this.pos.y.toFixed(3), 8)}`+
            `${f(this.pos.z.toFixed(3),8)}` +
            "  1.00  0.00              ";
        Atom.serial_atom++;
        if (Atom.serial_atom > 99999) {
            Atom.serial_atom = 1;
        }
        return res;
    }
}

export { Nucleotide, Atom, FROM_ANGSTROM_TO_OXDNA };
